from setuptools import setup, find_packages, Extension

# setup() parameters - https://packaging.python.org/guides/distributing-packages-using-setuptools/
setup(
    name='mynsc',
    version='0.1.0',
    author="name",
    description="nsc",
    install_requires=['invoke'],
    packages=find_packages(),

    license='GPL-2',

    ext_modules=[
        Extension(
            'mynsc',
            ['nsc.c'],
        ),
    ]
)
