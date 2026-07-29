from setuptools import setup, find_packages
import os

# Read README for long description
def read_readme():
    readme_file = os.path.join(os.path.dirname(__file__), 'README_PYTHON.md')
    if os.path.exists(readme_file):
        with open(readme_file, 'r', encoding='utf-8') as f:
            return f.read()
    return ''

setup(
    name='ceb-iv-fit',
    version='1.0.0',
    author='Scientific Computing Team',
    description='Cold Electron Bolometer IV Curve Fitting - Python Implementation',
    long_description=read_readme(),
    long_description_content_type='text/markdown',
    packages=find_packages(),
    python_requires='>=3.7',
    install_requires=[
        'numpy>=1.21.0',
        'scipy>=1.7.0',
        'lmfit>=1.0.0',
    ],
    extras_require={
        'dev': [
            'pytest>=6.0',
            'matplotlib>=3.0',
        ],
    },
    entry_points={
        'console_scripts': [
            'ceb-fit=main:main',
            'ceb-verify=verify_implementation:main',
            'ceb-convert=convert_config:main',
            'ceb-test=test_implementation:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    keywords='physics curve-fitting bolometers ceb iv-scientific',
    project_urls={
        'Documentation': 'https://github.com/yourusername/ceb_iv_fit',
        'Source': 'https://github.com/yourusername/ceb_iv_fit',
    },
)