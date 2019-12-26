from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='GNIRSLongSlit',
    version='0.3.0',
    packages = find_packages(),
    url='https://github.com/astrochun/GNIRSLongSlit',
    license='MIT License',
    author='Chun Ly',
    author_email='astro.chun@gmail.com',
    description='Python 2.7 codes to reduce Longslit data from Gemini-N/GNIRS',
    long_description = long_description,
    long_description_content_type = "text/markdown",
    install_requires = ['numpy', 'astropy']
)
