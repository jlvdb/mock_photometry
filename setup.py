import setuptools


with open("requirements.txt") as f:
    install_requires = [pkg.strip() for pkg in f.readlines() if pkg.strip()]

setuptools.setup(
    name="mock_photometry",
    version="1.0",
    author="Jan Luca van den Busch",
    description="Analytic prescription to simulate photometric noise in simulated galaxy catalogs.",
    long_description_content_type="text/markdown",
    url="https://github.com/jlvdb/mock_photometry",
    packages=setuptools.find_packages(),
    install_requires=install_requires)
