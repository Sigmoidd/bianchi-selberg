"""Exact finite images and adjacency operators, reusable without Arb."""
from .arithmetic import ResidueRing, ProjectiveMatrices
from .action import FiniteAction, generated_action, standard_generators
__all__ = ["ResidueRing", "ProjectiveMatrices", "FiniteAction", "generated_action", "standard_generators"]
