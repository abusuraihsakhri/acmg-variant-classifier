"""Security and defensive coding utilities.

Provides path traversal defense, Windows reserved device name neutralization,
CSV formula injection prevention, and strict numeric type and boundary validation.
"""

import math
import os
from pathlib import Path
from typing import Any, Union

WINDOWS_RESERVED_NAMES = frozenset({
    "CON", "PRN", "AUX", "NUL",
    *(f"COM{i}" for i in range(1, 10)),
    *(f"LPT{i}" for i in range(1, 10)),
})


def safe_resolve_path(path: Union[str, Path, os.PathLike], base_dir: Union[str, Path, os.PathLike, None] = None) -> Path:
    """Resolve a path safely, blocking null bytes, Windows reserved devices,
    and directory traversal outside base_dir.

    Args:
        path: User-supplied file or directory path.
        base_dir: Optional base directory to confine the path within.

    Returns:
        Canonical, resolved Path object.

    Raises:
        TypeError: If path is not a string or PathLike.
        ValueError: If path contains null bytes or references Windows reserved device names.
        PermissionError: If resolved path escapes base_dir.
    """
    if not isinstance(path, (str, Path, os.PathLike)):
        raise TypeError(f"Path must be string or PathLike, got {type(path).__name__}")

    path_str = str(path)
    if "\0" in path_str:
        raise ValueError("Path contains illegal null byte")

    # Check for Windows reserved device names in the stem and filename
    pure_path = Path(path_str)
    stem_upper = pure_path.stem.upper()
    name_upper = pure_path.name.upper()
    if stem_upper in WINDOWS_RESERVED_NAMES or name_upper in WINDOWS_RESERVED_NAMES:
        raise ValueError(f"Path references reserved system device name: {stem_upper}")

    resolved = Path(os.path.realpath(os.path.abspath(path_str)))

    if base_dir is not None:
        resolved_base = Path(os.path.realpath(os.path.abspath(str(base_dir))))
        try:
            common = Path(os.path.commonpath([str(resolved_base), str(resolved)]))
        except ValueError:
            raise PermissionError("Path traversal detected: cross-drive escape")
        if common != resolved_base:
            raise PermissionError(f"Path traversal detected: {path} escapes base directory {base_dir}")

    return resolved


def sanitize_csv_cell(val: Any) -> Any:
    """Neutralize CSV formula injection triggers (CWE-1236).

    Prepends a single quote to strings starting with '=', '+', '@', '\\t', '\\r',
    or '-' followed by non-numeric characters. Preserves legitimate numbers.
    """
    if val is None:
        return ""
    if isinstance(val, (int, float)) and not isinstance(val, bool):
        return val
    s = str(val)
    if not s:
        return ""
    if s[0] in ("=", "+", "@", "\t", "\r"):
        return f"'{s}"
    if s[0] == "-":
        try:
            float(s)
            return s  # Legitimate negative number like -4.5 or -1
        except ValueError:
            return f"'{s}"  # Dangerous text like -cmd|' /C calc'!A0
    return s


def validate_numeric_range(name: str, value: Any, low: float, high: float) -> None:
    """Validate that a numeric value is a real, finite number within bounds.

    Rejects booleans, NaN, Infinity, and values outside [low, high].

    Args:
        name: Parameter name for error reporting.
        value: The value to validate.
        low: Minimum allowable value (inclusive).
        high: Maximum allowable value (inclusive).

    Raises:
        TypeError: If value is a boolean or non-numeric.
        ValueError: If value is NaN, infinite, or outside [low, high].
    """
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be numeric, got {type(value).__name__}")
    if math.isnan(value):
        raise ValueError(f"{name} cannot be NaN")
    if math.isinf(value):
        raise ValueError(f"{name} cannot be infinite")
    if value < low or value > high:
        raise ValueError(f"{name} value {value} is outside valid range [{low}, {high}]")
