""" Install automol
"""
from distutils.core import setup


setup(name="automol",
      version="0.1.0",
      packages=["automol",
                "automol.inchi",
                "automol.graph",
                "automol.graph._dict",
                "automol.graph._inchi",
                "automol.graph._stereo",
                "automol.tests"],
      package_dir={'automol': "automol"},
      package_data={'automol': ["tests/data/*.txt"]},)
