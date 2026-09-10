from setuptools import setup

# Sphinx is needed only to BUILD THE DOCS, and importing it here made every install of
# this package depend on it -- `pip install -e . --no-deps` on a machine without sphinx
# failed before setup() was even reached. Register the command when it is available.
try:
    from sphinx.setup_command import BuildDoc
    _cmdclass = {'build_sphinx': BuildDoc}
except ImportError:
    BuildDoc = None
    _cmdclass = {}

name = 'spectral-tree-inference'
version = '0.1'
release = '0.1.0'
setup(
    name=name,
    version=version,
    description='Spectral methods for fitting Latent Tree models',
    url='https://github.com/NoahAmsel/spectral-tree-inference',
    author='Yariv Aizenbud, Noah Amsel, Ariel Jaffe, Amber Hu, Mamie Wang',
    author_email='yariv.aizenbud@yale.edu',
    license='GPLv3',
    packages=['spectraltree'],
    install_requires=[
        'dendropy',
        'numpy',
        'oct2py',
        'pandas',
        'pydantic>=2.0.0',
        'python-igraph',
        'scipy>=1.9.0',
        'scikit-learn',
        'seaborn',
        'sphinx',
        'sphinx-rtd-theme',
        'toytree',
        'tqdm'
    ],
    include_package_data=True,
    zip_safe=False,
    test_suite="tests",
    cmdclass=_cmdclass,
    command_options={
        'build_sphinx': {
            'project': ('setup.py', name),
            'version': ('setup.py', version),
            'release': ('setup.py', release),
            'source_dir': ('setup.py', 'docs'),
            'build_dir': ('setup.py', 'docs/_build')}})
