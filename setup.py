from setuptools import setup, Extension

setup(
    name='dnaio',
    version='0.1',
    author='',
    author_email='',
    url='',
    description='',
    long_description='',
    license='MIT',
    packages=['dnaio'],
    package_dir={'': 'src'},
    ext_modules=[Extension('dnaio._core', sources=['src/dnaio/_core.pyx'])],
    extras_require={
        'dev': ['Cython', 'pytest'],
    },
    install_requires=['xopen'],
    classifiers=[
            "Development Status :: 1 - Planning",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Natural Language :: English",
            "Programming Language :: Cython",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
