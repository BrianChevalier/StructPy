"""
StructPy is a pure Python structural analysis.
"""

__version__ = '0.2'

try:
	import cross_sections
	import structural_classes
except ImportError:
	import StructPy.cross_sections
	import StructPy.structural_classes