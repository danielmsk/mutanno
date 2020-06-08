pip uninstall -y mutanno
rm -rf build
rm -rf ./dist/*
python3 setup.py sdist bdist_wheel
pip install ./dist/mutanno-0.4.1-py3-none-any.whl
# ls -l /Users/pcaso/venv/lib/python3.7/site-packages/mutanno/
# ls -l /Users/pcaso/venv/lib/python3.7/site-packages/mutanno/web/
mutanno
