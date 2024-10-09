from setuptools import setup, find_packages


setup(
    name="BiSCAN_module",                # Nombre de tu módulo/paquete
    version="1.0",                   # Versión del paquete
    packages=find_packages(),        # Encuentra automáticamente todos los subdirectorios con __init__.py
    entry_points={
        'console_scripts': [
            'biscan=BiSCAN_module.biscan_main:main',  # Aquí se asume que en jaja.py tienes una función llamada main()
        ],
    },
)