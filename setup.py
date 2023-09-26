from setuptools import find_packages, setup, Extension

module = Extension("myspkmodule", ["spkmeansmodule.c", "spkmeans.c"])
modulekmeans = Extension("mykmeanssp", sources = ["kmeans.c"])

setup(
    name = "myspkmodule", 
    version = "0.1.0",
    author = "Itay Hirschel & Nicole Simoni",
    author_email = "example@email",
    description = "This is a package",
    install_requires = ['invoke'],
    packages = find_packages(),
    license = "GPL-2",
    ext_modules = [
        module,
        modulekmeans
    ]
)