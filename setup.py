""" Install automol
"""
from distutils.core import setup


setup(name="automol",
      version="0.1.1",
      packages=["automol",
                "automol._cnst",
                "automol.inchi",
                "automol.geom",
                "automol.zmatrix",
                "automol.graph",
                "automol.graph._dict",
                "automol.graph._inchi",
                "automol.graph._stereo",
                "automol.tests"],
      package_dir={'automol': "automol"},
      package_data={'automol': ["tests/data/*.txt"]},)
