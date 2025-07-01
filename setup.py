# setup.py
from setuptools import setup

setup(
    name='scBAMpler',
    version='0.1',
    py_modules=['scBAMpler'],
    entry_points={
        'console_scripts': ['scBAMpler = scBAMpler:main'],
    },
)
