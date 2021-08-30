from setuptools import setup, find_packages, Extension

# setup() parameters - https://packaging.python.org/guides/distributing-packages-using-setuptools/
setup(
    name='mykmeanssp',
    version='0.1.0',
    author="name",
    description="kmeanspp",
    install_requires=['invoke'],
    packages=find_packages(),

    license='GPL-2',

    ext_modules=[
        Extension(
            'mykmeanssp',
            ['kmeans.c'],
        ),
    ]
)
