#!/usr/bin/env python3
"""Check the compact magnetic-topology/QSL regression summary."""

import csv
import math
import os
import struct
import sys
import xml.etree.ElementTree as ET


TOL_UNIFORM = 1.0e-12
TOL_HYPERBOLIC_REL = 2.0e-6
TOL_FACE_LIMIT = 1.0e-12
TOL_SINGLE_MULTI = 1.0e-14
TOL_COMBINED = 1.0e-14
TOL_ARBITRARY_COMPARE = 1.0e-14
TOL_FIELDLINE_PRODUCTS = 1.0e-12


def parse_float(value):
    return float(value.strip())


def parse_int(value):
    return int(value.strip())


def load_rows(path):
    rows = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = (row["test_name"].strip(), row["subcase"].strip())
            rows[key] = {name: value.strip() for name, value in row.items()}
    return rows


def require(condition, message, failures):
    if condition:
        print(f"PASS {message}")
    else:
        print(f"FAIL {message}")
        failures.append(message)


def require_row(rows, test_name, subcase, failures):
    key = (test_name, subcase)
    require(key in rows, f"row exists: {test_name}/{subcase}", failures)
    return rows.get(key)


def load_header(path):
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        return next(reader)


def load_csv_dict_rows(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader)


def file_exists(path):
    return os.path.exists(path)


def load_vtu(path, failures):
    try:
        tree = ET.parse(path)
    except FileNotFoundError:
        require(False, f"VTU exists: {path}", failures)
        return None, None, set(), []
    except ET.ParseError as exc:
        require(False, f"VTU parses: {path}: {exc}", failures)
        return None, None, set(), []

    root = tree.getroot()
    require(root.tag == "VTKFile", f"VTU root tag: {path}", failures)
    require(root.get("type") == "UnstructuredGrid", f"VTU type=UnstructuredGrid: {path}", failures)
    piece = root.find("./UnstructuredGrid/Piece")
    require(piece is not None, f"VTU has Piece: {path}", failures)
    if piece is None:
        return root, None, set(), []

    names = {
        node.get("Name")
        for node in piece.findall(".//DataArray")
        if node.get("Name") is not None
    }
    types_node = piece.find("./Cells/DataArray[@Name='types']")
    cell_types = []
    if types_node is not None and types_node.text is not None:
        cell_types = [parse_int(value) for value in types_node.text.split()]
    return root, piece, names, cell_types


def load_vti_appended(path, failures):
    try:
        with open(path, "rb") as handle:
            content = handle.read()
    except FileNotFoundError:
        require(False, f"VTI exists: {path}", failures)
        return None, None, [], content if "content" in locals() else b"", -1

    appended_start = content.find(b"<AppendedData")
    require(appended_start >= 0, f"VTI has AppendedData: {path}", failures)
    if appended_start < 0:
        return None, None, [], content, -1

    header_xml = content[:appended_start].decode("utf-8", errors="strict") + "</VTKFile>"
    try:
        root = ET.fromstring(header_xml)
    except ET.ParseError as exc:
        require(False, f"VTI header parses: {path}: {exc}", failures)
        return None, None, [], content, -1

    require(root.tag == "VTKFile", f"VTI root tag: {path}", failures)
    require(root.get("type") == "ImageData", f"VTI type=ImageData: {path}", failures)
    image = root.find("./ImageData")
    require(image is not None, f"VTI has ImageData: {path}", failures)
    piece = root.find("./ImageData/Piece")
    require(piece is not None, f"VTI has Piece: {path}", failures)
    arrays = []
    if piece is not None:
        arrays = piece.findall("./PointData/DataArray")

    tag_end = content.find(b">", appended_start)
    payload_start = content.find(b"_", tag_end)
    require(payload_start >= 0, f"VTI has appended payload marker: {path}", failures)
    return root, image, arrays, content, payload_start + 1


def load_vti_payload(path, failures):
    root, image, arrays, content, payload_start = load_vti_appended(path, failures)
    values = {}
    array_info = {}
    if image is None or payload_start < 0:
        return root, image, array_info, values

    for node in arrays:
        name = node.get("Name")
        vtk_type = node.get("type")
        offset_text = node.get("offset")
        if name is None or vtk_type is None or offset_text is None:
            continue
        offset = parse_int(offset_text)
        byte_count = struct.unpack_from("<i", content, payload_start + offset)[0]
        start = payload_start + offset + 4
        if vtk_type == "Float64":
            nvalue = byte_count // 8
            values[name] = list(struct.unpack_from("<" + "d" * nvalue, content, start))
        elif vtk_type == "Int32":
            nvalue = byte_count // 4
            values[name] = list(struct.unpack_from("<" + "i" * nvalue, content, start))
        else:
            nvalue = byte_count
        array_info[name] = {
            "type": vtk_type,
            "offset": offset,
            "format": node.get("format"),
            "byte_count": byte_count,
            "nvalue": nvalue,
        }
    return root, image, array_info, values


def parse_extent(value):
    return [parse_int(part) for part in value.split()]


def parse_float_triplet(value):
    return [parse_float(part) for part in value.split()]


def load_vtu_point_data(path, failures):
    root, piece, names, cell_types = load_vtu(path, failures)
    data = {}
    if piece is None:
        return root, piece, names, cell_types, data

    point_data = piece.find("./PointData")
    if point_data is None:
        require(False, f"VTU has PointData: {path}", failures)
        return root, piece, names, cell_types, data

    for node in point_data.findall("./DataArray"):
        name = node.get("Name")
        if name is None:
            continue
        values = (node.text or "").split()
        if node.get("type") in ("Int32", "UInt8"):
            data[name] = [parse_int(value) for value in values]
        else:
            data[name] = [parse_float(value) for value in values]
    return root, piece, names, cell_types, data


def require_vtu_arrays(names, expected_names, label, failures):
    for name in expected_names:
        require(name in names, f"{label} has {name}", failures)


def check_uniform_axis_vtu(path, label, failures):
    required_names = (
        "length_total",
        "twist_total",
        "logQ",
        "connection_type_Q",
        "logQperp",
    )
    root, piece, names, cell_types, data = load_vtu_point_data(path, failures)
    if piece is None:
        return

    require(parse_int(piece.get("NumberOfPoints")) == 25, f"{label} VTU points=25", failures)
    require(parse_int(piece.get("NumberOfCells")) == 16, f"{label} VTU cells=16", failures)
    require(cell_types == [9] * 16, f"{label} VTU cell types are VTK_QUAD", failures)
    require_vtu_arrays(names, required_names, label, failures)

    forbidden = (
        "q_forward_" + "method" + "1",
        "q_backward_" + "method" + "1",
        "qperp_method2",
        "N2_qperp",
        "valid_Q",
        "status_Q",
        "face_pair_Q",
        "valid_Qperp",
        "status_Qperp",
        "logQ_" + "vis",
        "valid_Q_" + "vis",
        "q_" + "scott",
        "logq_" + "scott",
        "q_" + "hy" + "brid",
        "logQ_" + "hy" + "brid",
    )
    for old_name in forbidden:
        require(old_name not in names, f"{label} minimal VTU has no forbidden field {old_name}", failures)

    if "logQ" in data:
        log_values = data["logQ"]
        require(
            len(log_values) == 25 and max(abs(q - math.log10(2.0)) for q in log_values) < TOL_UNIFORM,
            f"{label} Scott q0 logQ=log10(2)",
            failures,
        )

    if "logQperp" in data:
        logqperp_values = data["logQperp"]
        require(
            len(logqperp_values) == 25 and max(abs(q - math.log10(2.0)) for q in logqperp_values) < TOL_UNIFORM,
            f"{label} logQperp=log10(2)",
            failures,
        )


def check_hyperbolic_axis_vtu(path, failures):
    expected = 2.0 * math.cosh(0.4)
    root, piece, names, cell_types, data = load_vtu_point_data(path, failures)
    if piece is None:
        return

    require(parse_int(piece.get("NumberOfPoints")) == 49, "xy hyperbolic VTU points=49", failures)
    require(parse_int(piece.get("NumberOfCells")) == 36, "xy hyperbolic VTU cells=36", failures)
    require(cell_types == [9] * 36, "xy hyperbolic VTU cell types are VTK_QUAD", failures)
    require("logQperp" in names, "xy hyperbolic VTU has logQperp", failures)
    require("logQ" in names, "xy hyperbolic VTU has logQ", failures)
    require("connection_type_Q" in names, "xy hyperbolic VTU has connection_type_Q", failures)
    require("qperp_method2" not in names, "xy hyperbolic minimal VTU has no qperp_method2", failures)
    require("q_forward_" + "method" + "1" not in names, "xy hyperbolic minimal VTU has no old forward-Q field", failures)
    require("logQ_" + "vis" not in names, "xy hyperbolic VTU has no old Q visual alias", failures)
    require("q_display" not in names, "xy hyperbolic VTU has no q_display", failures)

    if "logQperp" in data:
        center_index = 24
        logqperp_center = data["logQperp"][center_index]
        rel_err = abs(10.0 ** logqperp_center - expected) / expected
        require(rel_err < TOL_HYPERBOLIC_REL, "xy hyperbolic center logQperp rel_err", failures)
    if "logQ" in data:
        center_index = 24
        logq_center = data["logQ"][center_index]
        rel_err = abs(10.0 ** logq_center - expected) / expected
        require(rel_err < TOL_HYPERBOLIC_REL, "xy hyperbolic center Scott q0 logQ rel_err", failures)


def check_qperp_uniform(rows, failures):
    for plane in ("xy", "xz", "yz"):
        row = require_row(rows, "qperp_uniform", plane, failures)
        if row is None:
            continue
        require(parse_int(row["valid_count"]) == 25, f"Qperp uniform {plane} valid=25", failures)
        require(parse_int(row["invalid_count"]) == 0, f"Qperp uniform {plane} invalid=0", failures)
        require(parse_int(row["nonfinite_count"]) == 0, f"Qperp uniform {plane} nonfinite=0", failures)
        require(abs(parse_float(row["q_min"]) - 2.0) < TOL_UNIFORM, f"Qperp uniform {plane} q_min=2", failures)
        require(abs(parse_float(row["q_max"]) - 2.0) < TOL_UNIFORM, f"Qperp uniform {plane} q_max=2", failures)
        require(abs(parse_float(row["q_mean"]) - 2.0) < TOL_UNIFORM, f"Qperp uniform {plane} q_mean=2", failures)


def check_qperp_hyperbolic(rows, failures):
    row = require_row(rows, "qperp_hyperbolic_center", "center", failures)
    if row is None:
        return
    require(parse_int(row["valid_count"]) == 1, "Qperp hyperbolic center valid", failures)
    require(parse_float(row["rel_err"]) < TOL_HYPERBOLIC_REL, "Qperp hyperbolic center rel_err", failures)


def check_qperp_arbitrary(rows, failures):
    row = require_row(rows, "qperp_arbitrary", "eq_xy", failures)
    if row is not None:
        require(parse_int(row["valid_count"]) == 25, "arbitrary Qperp eq_xy valid=25", failures)
        require(parse_int(row["invalid_count"]) == 0, "arbitrary Qperp eq_xy invalid=0", failures)
        require(parse_int(row["nonfinite_count"]) == 0, "arbitrary Qperp eq_xy nonfinite=0", failures)
        require(abs(parse_float(row["q_min"]) - 2.0) < TOL_UNIFORM, "arbitrary Qperp eq_xy q_min=2", failures)
        require(abs(parse_float(row["q_max"]) - 2.0) < TOL_UNIFORM, "arbitrary Qperp eq_xy q_max=2", failures)
        require(abs(parse_float(row["q_mean"]) - 2.0) < TOL_UNIFORM, "arbitrary Qperp eq_xy q_mean=2", failures)
        require(parse_float(row["max_abs_diff"]) < TOL_ARBITRARY_COMPARE, "arbitrary Qperp eq_xy axis compare", failures)

    row = require_row(rows, "qperp_arbitrary", "tilted_uniform", failures)
    if row is not None:
        require(parse_int(row["valid_count"]) == 25, "arbitrary Qperp tilted valid=25", failures)
        require(parse_int(row["invalid_count"]) == 0, "arbitrary Qperp tilted invalid=0", failures)
        require(parse_int(row["nonfinite_count"]) == 0, "arbitrary Qperp tilted nonfinite=0", failures)
        require(abs(parse_float(row["q_min"]) - 2.0) < TOL_UNIFORM, "arbitrary Qperp tilted q_min=2", failures)
        require(abs(parse_float(row["q_max"]) - 2.0) < TOL_UNIFORM, "arbitrary Qperp tilted q_max=2", failures)
        require(abs(parse_float(row["q_mean"]) - 2.0) < TOL_UNIFORM, "arbitrary Qperp tilted q_mean=2", failures)

    row = require_row(rows, "qperp_arbitrary", "outside_seed", failures)
    if row is not None:
        require(parse_int(row["invalid_count"]) > 0, "arbitrary Qperp outside invalid_count>0", failures)
        require(parse_int(row["status_expected"]) == parse_int(row["status_actual"]), "arbitrary Qperp outside seed status", failures)


def check_face_limit(rows, failures):
    for face in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax"):
        row = require_row(rows, "endpoint_face_limit", face, failures)
        if row is None:
            continue
        require(parse_int(row["status_expected"]) == parse_int(row["status_actual"]), f"face-limit {face} status", failures)
        require(parse_float(row["abs_err"]) < TOL_FACE_LIMIT, f"face-limit {face} B error", failures)


def check_single_multi(rows, failures):
    row = require_row(rows, "qperp_single_multi", "nseed1", failures)
    if row is not None:
        require(parse_int(row["valid_count"]) == 1, "Qperp nseed=1 valid", failures)
        require(parse_int(row["status_expected"]) == parse_int(row["status_actual"]), "Qperp nseed=1 status", failures)
        require(parse_float(row["max_abs_diff"]) < TOL_SINGLE_MULTI, "Qperp single/multi max_abs_diff", failures)
    row = require_row(rows, "qperp_single_multi", "invalid_seed", failures)
    if row is not None:
        require(parse_int(row["valid_count"]) == 0, "Qperp invalid seed valid=0", failures)
        require(parse_int(row["invalid_count"]) == 1, "Qperp invalid seed invalid=1", failures)
        require(parse_int(row["status_expected"]) == parse_int(row["status_actual"]), "Qperp invalid seed status", failures)
        require(row["extra"] == "nan_scalars=1", "Qperp invalid seed NaN scalars", failures)


def check_fieldline_products_seeds(rows, failures):
    row = require_row(rows, "fieldline_products_seeds", "length_uniform", failures)
    if row is not None:
        require(parse_int(row["valid_count"]) == 3, "fieldline seeds length valid=3", failures)
        require(parse_int(row["invalid_count"]) == 1, "fieldline seeds length invalid=1", failures)
        require(parse_float(row["abs_err"]) < TOL_FIELDLINE_PRODUCTS, "fieldline seeds length error", failures)
        require(parse_float(row["max_abs_diff"]) < TOL_SINGLE_MULTI, "fieldline seeds length low-level diff", failures)
        require(parse_int(row["status_expected"]) == parse_int(row["status_actual"]), "fieldline seeds length outside status", failures)

    row = require_row(rows, "fieldline_products_seeds", "all_uniform", failures)
    if row is not None:
        require(parse_int(row["valid_count"]) == 3, "fieldline seeds all length valid=3", failures)
        require(parse_int(row["invalid_count"]) == 1, "fieldline seeds all length invalid=1", failures)
        require(abs(parse_float(row["q_min"]) - 2.0) < TOL_UNIFORM, "fieldline seeds all qperp q_min=2", failures)
        require(abs(parse_float(row["q_max"]) - 2.0) < TOL_UNIFORM, "fieldline seeds all qperp q_max=2", failures)
        require(abs(parse_float(row["q_mean"]) - 2.0) < TOL_UNIFORM, "fieldline seeds all qperp q_mean=2", failures)
        require(parse_float(row["abs_err"]) < TOL_FIELDLINE_PRODUCTS, "fieldline seeds all qperp error", failures)
        require(parse_float(row["max_abs_diff"]) < TOL_SINGLE_MULTI, "fieldline seeds all low-level diff", failures)
        require("twist_valid=3" in row["extra"], "fieldline seeds all twist valid=3", failures)
        require("qperp_valid=3" in row["extra"], "fieldline seeds all qperp valid=3", failures)
        require("twist_err= 0.00000E+00" in row["extra"], "fieldline seeds all twist error=0", failures)

    row = require_row(rows, "fieldline_products_seeds", "hyperbolic_center", failures)
    if row is not None:
        require(parse_int(row["valid_count"]) == 1, "fieldline seeds hyperbolic length valid=1", failures)
        require(abs(parse_float(row["q_mean"]) - parse_float(row["expected"])) < 1.0e-5, "fieldline seeds hyperbolic finite qperp", failures)
        require(parse_float(row["rel_err"]) < TOL_HYPERBOLIC_REL, "fieldline seeds hyperbolic rel_err", failures)
        require(parse_float(row["max_abs_diff"]) < TOL_SINGLE_MULTI, "fieldline seeds hyperbolic low-level diff", failures)

    length_header = load_header("reg_fieldline_products_seeds_length.csv")
    all_header = load_header("reg_fieldline_products_seeds_all.csv")
    require("twist_total" not in length_header, "fieldline seeds length-only header has no twist", failures)
    require("qperp" not in length_header, "fieldline seeds length-only header has no qperp", failures)
    require("twist_total" in all_header, "fieldline seeds all header has twist", failures)
    require("qperp" in all_header, "fieldline seeds all header has qperp", failures)


def check_fieldline_products_plane_arbitrary(rows, failures):
    row = require_row(rows, "fieldline_products_plane_arbitrary", "eq_xy_all", failures)
    if row is not None:
        require(parse_int(row["valid_count"]) == 25, "fieldline plane eq_xy length valid=25", failures)
        require(parse_int(row["invalid_count"]) == 0, "fieldline plane eq_xy length invalid=0", failures)
        require(abs(parse_float(row["q_min"]) - 2.0) < TOL_UNIFORM, "fieldline plane eq_xy q_min=2", failures)
        require(abs(parse_float(row["q_max"]) - 2.0) < TOL_UNIFORM, "fieldline plane eq_xy q_max=2", failures)
        require(abs(parse_float(row["q_mean"]) - 2.0) < TOL_UNIFORM, "fieldline plane eq_xy q_mean=2", failures)
        require(parse_float(row["abs_err"]) < TOL_FIELDLINE_PRODUCTS, "fieldline plane eq_xy qperp error", failures)
        require(parse_float(row["max_abs_diff"]) < TOL_SINGLE_MULTI, "fieldline plane eq_xy low-level diff", failures)
        require("twist_valid=25" in row["extra"], "fieldline plane eq_xy twist valid=25", failures)
        require("qperp_valid=25" in row["extra"], "fieldline plane eq_xy qperp valid=25", failures)
        require("twist_err= 0.00000E+00" in row["extra"], "fieldline plane eq_xy twist error=0", failures)

    row = require_row(rows, "fieldline_products_plane_arbitrary", "length_only_header", failures)
    if row is not None:
        require("header_has_twist=0" in row["extra"], "fieldline plane length-only summary no twist", failures)
        require("header_has_qperp=0" in row["extra"], "fieldline plane length-only summary no qperp", failures)

    row = require_row(rows, "fieldline_products_plane_arbitrary", "tilted_uniform", failures)
    if row is not None:
        require(parse_int(row["valid_count"]) == 25, "fieldline plane tilted length valid=25", failures)
        require(abs(parse_float(row["q_min"]) - 2.0) < TOL_UNIFORM, "fieldline plane tilted q_min=2", failures)
        require(abs(parse_float(row["q_max"]) - 2.0) < TOL_UNIFORM, "fieldline plane tilted q_max=2", failures)
        require(abs(parse_float(row["q_mean"]) - 2.0) < TOL_UNIFORM, "fieldline plane tilted q_mean=2", failures)
        require("twist_valid=25" in row["extra"], "fieldline plane tilted twist valid=25", failures)
        require("qperp_valid=25" in row["extra"], "fieldline plane tilted qperp valid=25", failures)
        require("twist_err= 0.00000E+00" in row["extra"], "fieldline plane tilted twist error=0", failures)

    row = require_row(rows, "fieldline_products_plane_arbitrary", "hyperbolic_center", failures)
    if row is not None:
        require(parse_int(row["valid_count"]) == 49, "fieldline plane hyperbolic length valid=49", failures)
        require(parse_float(row["rel_err"]) < TOL_HYPERBOLIC_REL, "fieldline plane hyperbolic rel_err", failures)
        require(parse_float(row["max_abs_diff"]) < TOL_SINGLE_MULTI, "fieldline plane hyperbolic low-level diff", failures)

    row = require_row(rows, "fieldline_products_plane_arbitrary", "outside_seed", failures)
    if row is not None:
        require(parse_int(row["invalid_count"]) > 0, "fieldline plane outside invalid length count", failures)
        require(parse_int(row["status_expected"]) == parse_int(row["status_actual"]), "fieldline plane outside status", failures)
        require("outside=" in row["extra"], "fieldline plane outside count present", failures)

    row = require_row(rows, "fieldline_products_plane_arbitrary", "invalid_basis", failures)
    if row is not None:
        require("header_only=1" in row["extra"], "fieldline plane invalid basis header-only", failures)

    length_header = load_header("reg_fieldline_products_plane_arbitrary_length.csv")
    all_header = load_header("reg_fieldline_products_plane_arbitrary_eq_xy_all.csv")
    require(length_header[:7] == ["i", "j", "s1", "s2", "seed_x", "seed_y", "seed_z"], "fieldline plane length header prefix", failures)
    require("twist_total" not in length_header, "fieldline plane length header has no twist", failures)
    require("qperp" not in length_header, "fieldline plane length header has no qperp", failures)
    require("twist_total" in all_header, "fieldline plane all header has twist", failures)
    require("qperp" in all_header, "fieldline plane all header has qperp", failures)


def check_vtu_smoke(failures):
    root, piece, names, cell_types = load_vtu("reg_fieldline_products_seeds_all.vtu", failures)
    if piece is not None:
        require(parse_int(piece.get("NumberOfPoints")) == 4, "VTU seed-set points=4", failures)
        require(parse_int(piece.get("NumberOfCells")) == 4, "VTU seed-set cells=4", failures)
        require(cell_types == [1, 1, 1, 1], "VTU seed-set cell types are VTK_VERTEX", failures)
        for name in ("length_total", "twist_total", "logQperp", "valid_Qperp", "status_Qperp"):
            require(name in names, f"VTU seed-set has {name}", failures)
        for name in ("qperp", "logqperp", "N2", "bfactor"):
            require(name not in names, f"VTU seed-set minimal has no forbidden field {name}", failures)

    root, piece, names, cell_types = load_vtu("reg_fieldline_products_plane_arbitrary_eq_xy_all.vtu", failures)
    if piece is not None:
        require(parse_int(piece.get("NumberOfPoints")) == 25, "VTU plane all points=25", failures)
        require(parse_int(piece.get("NumberOfCells")) == 16, "VTU plane all cells=16", failures)
        require(cell_types == [9] * 16, "VTU plane all cell types are VTK_QUAD", failures)
        for name in ("length_total", "twist_total", "logQperp", "valid_Qperp", "status_Qperp"):
            require(name in names, f"VTU plane all has {name}", failures)
        for name in ("qperp", "logqperp", "N2", "bfactor"):
            require(name not in names, f"VTU plane all minimal has no forbidden field {name}", failures)

    root, piece, names, cell_types = load_vtu("reg_fieldline_products_plane_arbitrary_length.vtu", failures)
    if piece is not None:
        require("length_total" in names, "VTU plane length-only has length_total", failures)
        require("twist_total" not in names, "VTU plane length-only has no twist_total", failures)
        require("logQperp" not in names, "VTU plane length-only has no logQperp", failures)

    root, piece, names, cell_types = load_vtu("reg_fieldline_products_plane_arbitrary_invalid_basis.vtu", failures)
    if piece is not None:
        require(parse_int(piece.get("NumberOfPoints")) == 0, "VTU invalid-basis points=0", failures)
        require(parse_int(piece.get("NumberOfCells")) == 0, "VTU invalid-basis cells=0", failures)

    check_uniform_axis_vtu("reg_qsl_plane_vtu_xy_uniform.vtu", "xy full", failures)
    check_uniform_axis_vtu("reg_qsl_plane_vtu_xz_uniform.vtu", "xz full", failures)
    check_uniform_axis_vtu("reg_qsl_plane_vtu_yz_uniform.vtu", "yz full", failures)
    check_hyperbolic_axis_vtu("reg_qsl_plane_vtu_xy_hyperbolic.vtu", failures)


def check_vti_smoke(rows, failures):
    row = require_row(rows, "vti_smoke", "pointdata", failures)
    if row is not None:
        require(parse_int(row["valid_count"]) == 24, "VTI smoke summary point count=24", failures)
        require(row["extra"] == "appended_binary", "VTI smoke summary appended_binary", failures)

    root, image, arrays, content, payload_start = load_vti_appended(
        "reg_vti_smoke_pointdata.vti", failures
    )
    if image is None or payload_start < 0:
        return

    require(parse_extent(image.get("WholeExtent")) == [0, 3, 0, 2, 0, 1], "VTI WholeExtent", failures)
    require(
        all(abs(a - b) < 1.0e-12 for a, b in zip(parse_float_triplet(image.get("Origin")), [-0.3, -0.2, -0.1])),
        "VTI Origin",
        failures,
    )
    require(
        all(abs(a - b) < 1.0e-12 for a, b in zip(parse_float_triplet(image.get("Spacing")), [0.2, 0.2, 0.2])),
        "VTI Spacing",
        failures,
    )

    array_by_name = {node.get("Name"): node for node in arrays}
    for name in ("length_total", "qperp", "status"):
        require(name in array_by_name, f"VTI has {name}", failures)
        if name in array_by_name:
            require(array_by_name[name].get("format") == "appended", f"VTI {name} is appended", failures)
            require(array_by_name[name].get("offset") is not None, f"VTI {name} has offset", failures)

    if not all(name in array_by_name for name in ("length_total", "qperp", "status")):
        return
    offsets = [parse_int(array_by_name[name].get("offset")) for name in ("length_total", "qperp", "status")]
    require(offsets == sorted(offsets), "VTI offsets monotonically increase", failures)
    require(offsets == [0, 196, 392], "VTI offsets match expected byte layout", failures)

    expected_byte_counts = {
        "length_total": 24 * 8,
        "qperp": 24 * 8,
        "status": 24 * 4,
    }
    for name in ("length_total", "qperp", "status"):
        offset = parse_int(array_by_name[name].get("offset"))
        byte_count = struct.unpack_from("<i", content, payload_start + offset)[0]
        require(byte_count == expected_byte_counts[name], f"VTI {name} byte count", failures)


def require_vti_payload_arrays(array_info, expected_names, forbidden_names, label, failures):
    for name in expected_names:
        require(name in array_info, f"{label} VTI has {name}", failures)
        if name in array_info:
            require(array_info[name]["format"] == "appended", f"{label} VTI {name} is appended", failures)
    for name in forbidden_names:
        require(name not in array_info, f"{label} VTI has no forbidden field {name}", failures)


def require_vti_byte_counts(array_info, npoint, label, failures):
    offsets = []
    for name, info in array_info.items():
        offsets.append(info["offset"])
        expected = npoint * (8 if info["type"] == "Float64" else 4)
        require(info["byte_count"] == expected, f"{label} VTI {name} byte count", failures)
    require(offsets == sorted(offsets), f"{label} VTI offsets monotonically increase", failures)


def max_abs_diff(left, right):
    diff = 0.0
    for a, b in zip(left, right):
        if math.isnan(a) and math.isnan(b):
            continue
        diff = max(diff, abs(a - b))
    return diff


def max_int_diff(left, right):
    return max((abs(a - b) for a, b in zip(left, right)), default=0)


def check_volume_vti(rows, failures):
    row = require_row(rows, "volume_vti", "length_uniform", failures)
    root, image, array_info, values = load_vti_payload("reg_volume_vti_length_uniform.vti", failures)
    if image is not None:
        require(parse_extent(image.get("WholeExtent")) == [0, 3, 0, 2, 0, 1], "volume length VTI WholeExtent", failures)
        require_vti_payload_arrays(
            array_info,
            ("length_total",),
            ("twist_total", "qperp"),
            "volume length",
            failures,
        )
        require_vti_byte_counts(array_info, 24, "volume length", failures)
        if "length_total" in values:
            require(max(abs(value - 2.0) for value in values["length_total"]) < TOL_FIELDLINE_PRODUCTS, "volume length VTI length_total=2", failures)
    if row is not None:
        require(parse_int(row["valid_count"]) == 24, "volume length summary valid=24", failures)
        require(parse_int(row["invalid_count"]) == 0, "volume length summary invalid=0", failures)
        require(parse_float(row["abs_err"]) < TOL_FIELDLINE_PRODUCTS, "volume length summary error", failures)

    row = require_row(rows, "volume_vti", "all_uniform", failures)
    root, image, array_info, values = load_vti_payload("reg_volume_vti_all_uniform.vti", failures)
    if image is not None:
        require_vti_payload_arrays(
            array_info,
            (
                "length_total",
                "twist_total",
                "logQperp",
                "valid_Qperp",
                "status_Qperp",
            ),
            ("qperp", "logqperp", "N2", "bfactor", "face_forward_qperp"),
            "volume all",
            failures,
        )
        require_vti_byte_counts(array_info, 24, "volume all", failures)
        if "length_total" in values:
            require(max(abs(value - 2.0) for value in values["length_total"]) < TOL_FIELDLINE_PRODUCTS, "volume all VTI length_total=2", failures)
        if "twist_total" in values:
            require(max(abs(value) for value in values["twist_total"]) < TOL_FIELDLINE_PRODUCTS, "volume all VTI twist_total=0", failures)
        if "logQperp" in values:
            require(max(abs(value - math.log10(2.0)) for value in values["logQperp"]) < TOL_UNIFORM, "volume all VTI logQperp=log10(2)", failures)
        if "valid_Qperp" in values:
            require(sum(values["valid_Qperp"]) == 24, "volume all VTI valid_Qperp=24", failures)
    if row is not None:
        require(parse_int(row["valid_count"]) == 24, "volume all summary length valid=24", failures)
        require(abs(parse_float(row["q_min"]) - 2.0) < TOL_UNIFORM, "volume all summary q_min=2", failures)
        require(abs(parse_float(row["q_max"]) - 2.0) < TOL_UNIFORM, "volume all summary q_max=2", failures)
        require(parse_float(row["abs_err"]) < TOL_FIELDLINE_PRODUCTS, "volume all summary qperp error", failures)
        require(parse_float(row["max_abs_diff"]) < TOL_FIELDLINE_PRODUCTS, "volume all summary twist error", failures)
        require("qperp_valid=24" in row["extra"], "volume all summary qperp valid=24", failures)

    row = require_row(rows, "volume_vti", "slab_equivalence", failures)
    _, _, chunk_all_info, chunk_all_values = load_vti_payload(
        "reg_volume_vti_all_uniform_chunk_all.vti", failures
    )
    _, _, chunk1_info, chunk1_values = load_vti_payload(
        "reg_volume_vti_all_uniform_chunk1.vti", failures
    )
    if chunk_all_info and chunk1_info:
        require_vti_payload_arrays(
            chunk_all_info,
            ("length_total", "twist_total", "logQperp", "valid_Qperp", "status_Qperp"),
            ("qperp", "N2", "bfactor"),
            "volume slab chunk_all",
            failures,
        )
        require_vti_payload_arrays(
            chunk1_info,
            ("length_total", "twist_total", "logQperp", "valid_Qperp", "status_Qperp"),
            ("qperp", "N2", "bfactor"),
            "volume slab chunk1",
            failures,
        )
        require_vti_byte_counts(chunk_all_info, 24, "volume slab chunk_all", failures)
        require_vti_byte_counts(chunk1_info, 24, "volume slab chunk1", failures)
        for name in ("length_total", "twist_total", "logQperp"):
            if name in chunk_all_values and name in chunk1_values:
                require(
                    max_abs_diff(chunk_all_values[name], chunk1_values[name]) < TOL_FIELDLINE_PRODUCTS,
                    f"volume slab {name} chunk_all vs chunk1",
                    failures,
                )
        for name in ("valid_Qperp", "status_Qperp"):
            if name in chunk_all_values and name in chunk1_values:
                require(
                    max_int_diff(chunk_all_values[name], chunk1_values[name]) == 0,
                    f"volume slab {name} chunk_all vs chunk1",
                    failures,
                )
    if row is not None:
        require(row["extra"] == "chunk_all_vs_chunk1", "volume slab summary extra", failures)

    row = require_row(rows, "volume_vti", "hyperbolic_center", failures)
    root, image, array_info, values = load_vti_payload("reg_volume_vti_hyperbolic_center.vti", failures)
    if image is not None:
        require_vti_payload_arrays(array_info, ("logQperp", "valid_Qperp", "status_Qperp"), ("qperp",), "volume hyperbolic", failures)
        require_vti_byte_counts(array_info, 27, "volume hyperbolic", failures)
        if "logQperp" in values:
            expected = 2.0 * math.cosh(0.4)
            rel_err = abs(10.0 ** values["logQperp"][13] - expected) / expected
            require(rel_err < TOL_HYPERBOLIC_REL, "volume hyperbolic VTI center qperp rel_err", failures)
    if row is not None:
        require(parse_float(row["rel_err"]) < TOL_HYPERBOLIC_REL, "volume hyperbolic summary rel_err", failures)

    row = require_row(rows, "volume_vti", "hyperbolic_center_slabbed", failures)
    root, image, array_info, values = load_vti_payload("reg_volume_vti_hyperbolic_center_chunk1.vti", failures)
    if image is not None:
        require_vti_payload_arrays(array_info, ("logQperp", "valid_Qperp", "status_Qperp"), ("qperp",), "volume hyperbolic slabbed", failures)
        require_vti_byte_counts(array_info, 27, "volume hyperbolic slabbed", failures)
        if "logQperp" in values:
            expected = 2.0 * math.cosh(0.4)
            rel_err = abs(10.0 ** values["logQperp"][13] - expected) / expected
            require(rel_err < TOL_HYPERBOLIC_REL, "volume hyperbolic slabbed VTI center qperp rel_err", failures)
    if row is not None:
        require(parse_float(row["rel_err"]) < TOL_HYPERBOLIC_REL, "volume hyperbolic slabbed summary rel_err", failures)

    row = require_row(rows, "volume_vti", "outside", failures)
    root, image, array_info, values = load_vti_payload("reg_volume_vti_outside.vti", failures)
    if image is not None:
        require_vti_payload_arrays(array_info, ("logQperp", "valid_Qperp", "status_Qperp"), ("qperp",), "volume outside", failures)
        require_vti_byte_counts(array_info, 18, "volume outside", failures)
        if "status_Qperp" in values:
            require(sum(1 for value in values["status_Qperp"] if value == 4) == 6, "volume outside VTI seed_outside count=6", failures)
    if row is not None:
        require(parse_int(row["status_expected"]) == parse_int(row["status_actual"]), "volume outside summary status", failures)
        require("outside=6" in row["extra"], "volume outside summary count=6", failures)

    row = require_row(rows, "volume_vti", "outside_slabbed", failures)
    root, image, array_info, values = load_vti_payload("reg_volume_vti_outside_chunk1.vti", failures)
    if image is not None:
        require_vti_payload_arrays(array_info, ("logQperp", "valid_Qperp", "status_Qperp"), ("qperp",), "volume outside slabbed", failures)
        require_vti_byte_counts(array_info, 18, "volume outside slabbed", failures)
        if "status_Qperp" in values:
            require(sum(1 for value in values["status_Qperp"] if value == 4) == 6, "volume outside slabbed VTI seed_outside count=6", failures)
    if row is not None:
        require(parse_int(row["status_expected"]) == parse_int(row["status_actual"]), "volume outside slabbed summary status", failures)
        require("outside=6" in row["extra"], "volume outside slabbed summary count=6", failures)


def check_namelist_runner(rows, failures):
    row = require_row(rows, "namelist_runner", "axis_plane_full_vtu", failures)
    root, piece, names, cell_types, data = load_vtu_point_data(
        "reg_namelist_qsl_xy_full.vtu", failures
    )
    if piece is not None:
        require(parse_int(piece.get("NumberOfPoints")) == 25, "namelist axis VTU points=25", failures)
        require(parse_int(piece.get("NumberOfCells")) == 16, "namelist axis VTU cells=16", failures)
        require(cell_types == [9] * 16, "namelist axis VTU cell types are VTK_QUAD", failures)
        require_vtu_arrays(
            names,
            (
                "length_total",
                "twist_total",
                "logQ",
                "connection_type_Q",
                "logQperp",
            ),
            "namelist axis",
            failures,
        )
        require("q_forward_" + "method" + "1" not in names, "namelist axis minimal has no old forward-Q field", failures)
        require("qperp_method2" not in names, "namelist axis minimal has no qperp_method2", failures)
        for name in ("valid_Q", "status_Q", "face_pair_Q", "valid_Qperp", "status_Qperp"):
            require(name not in names, f"namelist axis minimal has no {name}", failures)
        require("q_display" not in names, "namelist axis VTU has no q_display", failures)
    if row is not None:
        require(row["extra"] == "reg_namelist_qsl_xy_full.vtu", "namelist axis summary output", failures)

    row = require_row(rows, "namelist_runner", "axis_plane_full_vtu_full_detail", failures)
    root, piece, names, cell_types, data = load_vtu_point_data(
        "reg_namelist_qsl_xy_full_detail.vtu", failures
    )
    if piece is not None:
        require(parse_int(piece.get("NumberOfPoints")) == 25, "namelist full-detail axis VTU points=25", failures)
        require(parse_int(piece.get("NumberOfCells")) == 16, "namelist full-detail axis VTU cells=16", failures)
        require(cell_types == [9] * 16, "namelist full-detail axis VTU cell types are VTK_QUAD", failures)
        require_vtu_arrays(
            names,
            (
                "q",
                "twist_total",
                "logQ",
                "valid_Q",
                "status_Q",
                "face_pair_Q",
                "connection_type_Q",
                "qperp_method2",
                "logqperp_method2",
                "valid_qperp_method2",
            ),
            "namelist full-detail axis",
            failures,
        )
        forbidden = (
            "q_" + "method" + "1",
            "logq_" + "method" + "1",
            "q_forward_" + "method" + "1",
            "logq_forward_" + "method" + "1",
            "suspect_" + "method" + "1",
            "logQ_" + "vis",
            "valid_Q_" + "vis",
            "q_" + "scott",
            "logq_" + "scott",
            "q_" + "hy" + "brid",
            "logq_" + "hy" + "brid",
        )
        for name in forbidden:
            require(name not in names, f"namelist full-detail axis has no forbidden field {name}", failures)
        require("q_display" not in names, "namelist full-detail axis VTU has no q_display", failures)
    if row is not None:
        require(row["extra"] == "reg_namelist_qsl_xy_full_detail.vtu", "namelist full-detail axis summary output", failures)

    row = require_row(rows, "namelist_runner", "axis_plane_full_vtu_prefix", failures)
    root, piece, names, cell_types, data = load_vtu_point_data(
        "reg_namelist_prefix_xy_full.vtu", failures
    )
    if piece is not None:
        require(parse_int(piece.get("NumberOfPoints")) == 25, "namelist prefix axis VTU points=25", failures)
        require(parse_int(piece.get("NumberOfCells")) == 16, "namelist prefix axis VTU cells=16", failures)
        require(cell_types == [9] * 16, "namelist prefix axis VTU cell types are VTK_QUAD", failures)
        require_vtu_arrays(
            names,
            (
                "length_total",
                "twist_total",
                "logQ",
                "connection_type_Q",
                "logQperp",
            ),
            "namelist prefix axis",
            failures,
        )
        require("q_forward_" + "method" + "1" not in names, "namelist prefix axis minimal has no old forward-Q field", failures)
        require("qperp_method2" not in names, "namelist prefix axis minimal has no qperp_method2", failures)
        for name in ("valid_Q", "status_Q", "face_pair_Q", "valid_Qperp", "status_Qperp"):
            require(name not in names, f"namelist prefix axis minimal has no {name}", failures)
        require("q_display" not in names, "namelist prefix axis VTU has no q_display", failures)
    if row is not None:
        require(row["extra"] == "reg_namelist_prefix_xy_full.vtu", "namelist prefix axis summary output", failures)

    row = require_row(rows, "namelist_runner", "volume_vti", failures)
    root, image, array_info, values = load_vti_payload("reg_namelist_volume.vti", failures)
    if image is not None:
        require(parse_extent(image.get("WholeExtent")) == [0, 3, 0, 2, 0, 1], "namelist volume VTI WholeExtent", failures)
        require_vti_payload_arrays(
            array_info,
            ("length_total", "twist_total", "logQperp", "valid_Qperp", "status_Qperp"),
            ("q_forward_" + "method" + "1", "q_backward_" + "method" + "1"),
            "namelist volume",
            failures,
        )
        require_vti_byte_counts(array_info, 24, "namelist volume", failures)
    if row is not None:
        require(row["extra"] == "reg_namelist_volume.vti", "namelist volume summary output", failures)

    row = require_row(rows, "namelist_runner", "arbitrary_plane", failures)
    root, piece, names, cell_types, data = load_vtu_point_data(
        "reg_namelist_arb_arbitrary_plane_products.vtu", failures
    )
    if piece is not None:
        require(parse_int(piece.get("NumberOfPoints")) == 25, "namelist arbitrary VTU points=25", failures)
        require(parse_int(piece.get("NumberOfCells")) == 16, "namelist arbitrary VTU cells=16", failures)
        require(cell_types == [9] * 16, "namelist arbitrary VTU cell types are VTK_QUAD", failures)
        require_vtu_arrays(
            names,
            ("length_total", "twist_total", "logQperp", "valid_Qperp", "status_Qperp"),
            "namelist arbitrary",
            failures,
        )
        require("q_forward_" + "method" + "1" not in names, "namelist arbitrary has no old forward-Q field", failures)
        require("q_backward_" + "method" + "1" not in names, "namelist arbitrary has no old backward-Q field", failures)
        require("q_display" not in names, "namelist arbitrary has no q_display", failures)
    try:
        with open("reg_namelist_arb_arbitrary_plane_products.csv", "rb"):
            arb_csv_exists = True
    except FileNotFoundError:
        arb_csv_exists = False
    require(not arb_csv_exists, "namelist arbitrary default creates no CSV", failures)
    if row is not None:
        require(parse_int(row["status_actual"]) == 0, "namelist arbitrary summary CSV absent", failures)
        require("csv_absent=1" in row["extra"], "namelist arbitrary summary extra", failures)

    row = require_row(rows, "namelist_runner", "arbitrary_plane_csv", failures)
    root, piece, names, cell_types, data = load_vtu_point_data(
        "reg_namelist_arb_csv_arbitrary_plane_products.vtu", failures
    )
    if piece is not None:
        require(parse_int(piece.get("NumberOfPoints")) == 25, "namelist arbitrary CSV VTU points=25", failures)
        require(parse_int(piece.get("NumberOfCells")) == 16, "namelist arbitrary CSV VTU cells=16", failures)
        require(cell_types == [9] * 16, "namelist arbitrary CSV VTU cell types are VTK_QUAD", failures)
        require_vtu_arrays(
            names,
            ("length_total", "twist_total", "logQperp", "valid_Qperp", "status_Qperp"),
            "namelist arbitrary CSV",
            failures,
        )
        require("q_forward_" + "method" + "1" not in names, "namelist arbitrary CSV has no old forward-Q field", failures)
        require("q_backward_" + "method" + "1" not in names, "namelist arbitrary CSV has no old backward-Q field", failures)
        require("q_display" not in names, "namelist arbitrary CSV has no q_display", failures)
    try:
        header = load_header("reg_namelist_arb_csv_arbitrary_plane_products.csv")
        arb_csv_present = True
    except FileNotFoundError:
        header = []
        arb_csv_present = False
    require(arb_csv_present, "namelist arbitrary optional CSV exists", failures)
    if arb_csv_present:
        require(header[:7] == ["i", "j", "s1", "s2", "seed_x", "seed_y", "seed_z"], "namelist arbitrary CSV header prefix", failures)
        require("twist_total" in header, "namelist arbitrary CSV has twist_total", failures)
        require("qperp" in header, "namelist arbitrary CSV has qperp", failures)
    if row is not None:
        require(parse_int(row["status_actual"]) == 1, "namelist arbitrary CSV summary present", failures)
        require("csv_present=1" in row["extra"], "namelist arbitrary CSV summary extra", failures)

    row = require_row(rows, "namelist_runner", "seed_products", failures)
    seed_rows = load_csv_dict_rows("reg_namelist_seed_seed_products.csv")
    seed_header = list(seed_rows[0].keys()) if seed_rows else []
    require(seed_header[:4] == ["seed_id", "seed_x", "seed_y", "seed_z"], "namelist seed CSV header prefix", failures)
    require(len(seed_rows) == 4, "namelist seed CSV row count=4", failures)
    require(
        [parse_int(seed_row["seed_id"]) for seed_row in seed_rows] == [1, 2, 3, 4],
        "namelist seed CSV seed_id order",
        failures,
    )
    require("twist_total" in seed_header, "namelist seed CSV has twist_total", failures)
    require("qperp" in seed_header, "namelist seed CSV has qperp", failures)
    require("valid_qperp" in seed_header, "namelist seed CSV has valid_qperp", failures)
    require("status_qperp" in seed_header, "namelist seed CSV has status_qperp", failures)
    require(
        any(parse_int(seed_row["status_qperp"]) == 4 for seed_row in seed_rows),
        "namelist seed CSV has outside status=4",
        failures,
    )
    root, piece, names, cell_types, data = load_vtu_point_data(
        "reg_namelist_seed_seed_products.vtu", failures
    )
    if piece is not None:
        require(parse_int(piece.get("NumberOfPoints")) == 4, "namelist seed VTU points=4", failures)
        require(parse_int(piece.get("NumberOfCells")) == 4, "namelist seed VTU cells=4", failures)
        require(cell_types == [1] * 4, "namelist seed VTU cell types are VTK_VERTEX", failures)
        require_vtu_arrays(
            names,
            ("length_total", "twist_total", "logQperp", "valid_Qperp", "status_Qperp"),
            "namelist seed",
            failures,
        )
        require("q_forward_" + "method" + "1" not in names, "namelist seed has no old forward-Q field", failures)
        require("q_backward_" + "method" + "1" not in names, "namelist seed has no old backward-Q field", failures)
        require("q_display" not in names, "namelist seed has no q_display", failures)
    if row is not None:
        require(parse_int(row["status_actual"]) == 1, "namelist seed summary outputs present", failures)
        require("vtu_present=1" in row["extra"], "namelist seed summary extra", failures)

    row = require_row(rows, "namelist_runner", "axis_plane_csv_default", failures)
    axis_default_length = "reg_axis_csv_default_xy_length.csv"
    axis_default_optional = (
        "reg_axis_csv_default_xy_twist.csv",
        "reg_axis_csv_default_xy_mapping.csv",
        "reg_axis_csv_default_xy_qperp_method2.csv",
    )
    require(file_exists(axis_default_length), "namelist axis CSV default length exists", failures)
    for path in axis_default_optional:
        require(not file_exists(path), f"namelist axis CSV default optional absent: {path}", failures)
    if file_exists(axis_default_length):
        length_rows = load_csv_dict_rows(axis_default_length)
        length_header = list(length_rows[0].keys()) if length_rows else load_header(axis_default_length)
        require(len(length_rows) == 25, "namelist axis CSV default length row count=25", failures)
        require("length_total" in length_header, "namelist axis CSV default has length_total", failures)
        require("q_display" not in length_header, "namelist axis CSV default has no q_display", failures)
    if row is not None:
        require(parse_int(row["status_actual"]) == 1, "namelist axis CSV default summary ok", failures)
        require("optional_absent=1" in row["extra"], "namelist axis CSV default summary extra", failures)

    row = require_row(rows, "namelist_runner", "axis_plane_csv_all", failures)
    axis_all_files = (
        "reg_axis_csv_all_xy_length.csv",
        "reg_axis_csv_all_xy_twist.csv",
        "reg_axis_csv_all_xy_qperp_method2.csv",
    )
    for path in axis_all_files:
        require(file_exists(path), f"namelist axis CSV all file exists: {path}", failures)
    if all(file_exists(path) for path in axis_all_files):
        length_header = load_header("reg_axis_csv_all_xy_length.csv")
        twist_header = load_header("reg_axis_csv_all_xy_twist.csv")
        qperp_header = load_header("reg_axis_csv_all_xy_qperp_method2.csv")
        length_rows = load_csv_dict_rows("reg_axis_csv_all_xy_length.csv")
        require(len(length_rows) == 25, "namelist axis CSV all length row count=25", failures)
        require("length_total" in length_header, "namelist axis CSV all has length_total", failures)
        require("twist_total" in twist_header, "namelist axis CSV all has twist_total", failures)
        require("qperp" in qperp_header, "namelist axis CSV all Qperp has qperp", failures)
        require("valid" in qperp_header, "namelist axis CSV all Qperp has valid", failures)
        require("status" in qperp_header, "namelist axis CSV all Qperp has status", failures)
        for label, header in (
            ("length", length_header),
            ("twist", twist_header),
            ("Qperp", qperp_header),
        ):
            require("q_display" not in header, f"namelist axis CSV all {label} has no q_display", failures)
    if row is not None:
        require(parse_int(row["status_actual"]) == 1, "namelist axis CSV all summary ok", failures)
        require("all_present=1" in row["extra"], "namelist axis CSV all summary extra", failures)

    row = require_row(rows, "namelist_runner", "disabled", failures)
    if row is not None:
        require(parse_int(row["status_actual"]) == 0, "namelist disabled output absent summary", failures)
        require(row["extra"] == "output_absent=1", "namelist disabled extra", failures)
    try:
        with open("reg_namelist_disabled.vti", "rb"):
            disabled_exists = True
    except FileNotFoundError:
        disabled_exists = False
    require(not disabled_exists, "namelist disabled creates no output", failures)


def check_combined(rows, failures):
    for subcase in ("length", "twist", "mapping"):
        row = require_row(rows, "topology_combined_xy", subcase, failures)
        if row is None:
            continue
        require(parse_float(row["max_abs_diff"]) < TOL_COMBINED, f"combined xy {subcase} exact compare", failures)


def main(argv):
    if len(argv) != 2:
        print("usage: check_topology_qsl.py qsl_regression_summary.csv", file=sys.stderr)
        return 2

    rows = load_rows(argv[1])
    failures = []
    check_qperp_uniform(rows, failures)
    check_qperp_arbitrary(rows, failures)
    check_qperp_hyperbolic(rows, failures)
    check_face_limit(rows, failures)
    check_single_multi(rows, failures)
    check_fieldline_products_seeds(rows, failures)
    check_fieldline_products_plane_arbitrary(rows, failures)
    check_vtu_smoke(failures)
    check_vti_smoke(rows, failures)
    check_volume_vti(rows, failures)
    check_namelist_runner(rows, failures)
    check_combined(rows, failures)

    if failures:
        print(f"FAIL {len(failures)} check(s) failed")
        return 1
    print("PASS all magnetic topology QSL regression checks")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
