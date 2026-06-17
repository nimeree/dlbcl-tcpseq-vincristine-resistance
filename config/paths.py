"""Portable path helpers for thesis analysis scripts."""
from pathlib import Path
import os


def _root(env_var: str, default: str) -> Path:
    return Path(os.environ.get(env_var, default)).expanduser()


def analysis_path(*parts: str) -> Path:
    return _root(THESIS_ANALYSIS_DIR, data/analysis).joinpath(*parts)


def input_path(*parts: str) -> Path:
    return _root(THESIS_INPUT_DIR, data/input).joinpath(*parts)


def project_resource_path(*parts: str) -> Path:
    return _root(THESIS_RESOURCE_DIR, data/external).joinpath(*parts)


def external_path(*parts: str) -> Path:
    return _root(THESIS_EXTERNAL_ROOT, data/external).joinpath(*parts)
