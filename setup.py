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
    ext_modules=[Extension('_seqio', sources=['src/dnaio/_seqio.pyx'])],
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
