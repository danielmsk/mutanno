pip uninstall -y mutanno
rm -rf build
rm -rf ./dist/*
python3 setup.py sdist bdist_wheel
pip install ./dist/mutanno-0.3.14-py3-none-any.whl
mutanno
