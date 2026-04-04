from setuptools import setup, find_packages

setup(
    name='concretedesignpy',
    version='0.4.0',
    description='Open-source reinforced concrete design library with web interface (ACI 318-19, NSCP 2015, ACI 440R-17)',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    author='Albert Pamonag',
    author_email='albert@apeconsultancy.net',
    url='https://github.com/albertp16/concretedesignpy',
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'concretedesignpy': [
            'webapp/templates/**/*.html',
            'webapp/static/**/*',
        ],
    },
    install_requires=[
        'matplotlib',
        'numpy',
        'flask>=2.0',
        'gunicorn',
        'openpyxl',
    ],
    extras_require={
        'dev': ['pytest'],
    },
    entry_points={
        'console_scripts': [
            'concretedesignpy-web=concretedesignpy.webapp.app:main',
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    python_requires='>=3.8',
)