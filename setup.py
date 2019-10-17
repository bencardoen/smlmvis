from setuptools import setup
import versioneer

requirements = [
    'numpy', 'vtk', 'jupyter', 'pandas', 'requests', 'seaborn', 'scipy'
]

setup(
    name='smlmvis',
    license='MIT',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Superresolution visualization of 3D protein localization data from a range of microscopes",
    author="Ben Cardoen",
    author_email='bcardoen@sfu.ca',
    url='https://github.com/bencardoen/smlmvis',
    download_url='https://github.com/bencardoen/smlmvis/archive/v0.0.8.tar.gz',
    packages=['smlmvis'],
    install_requires=requirements,
    keywords='smlmvis',
    include_package_data=True,
    package_data={'': ['versioneer.py']},
    classifiers=[
        'Programming Language :: Python :: 3.7',
    ]
)
