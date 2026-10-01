"""Group interfaces, loaded lazily so exact arithmetic stays independent."""
from importlib import import_module

_EXPORTS = {"EllipticClass": "data", "GroupData": "data", "get_group": "data",
            "GroupKey": "identity", "LevelIdeal": "identity", "GroupRegistry": "registry",
            "InventoryBackend": "registry", "default_registry": "registry"}
__all__ = list(_EXPORTS)


def __getattr__(name):
    if name not in _EXPORTS:
        raise AttributeError(name)
    value = getattr(import_module(f"groups.{_EXPORTS[name]}"), name)
    globals()[name] = value
    return value
