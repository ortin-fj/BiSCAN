from setuptools import setup, find_packages


setup(
    name="BiSCAN_module",            
    version="1.0",                  
    packages=find_packages(),      
    entry_points={
        'console_scripts': [
            'biscan=BiSCAN_module.biscan_main:main', 
        ],
    },
)
