#https://qiita.com/shonansurvivors/items/0fbcbfde129f2d26301c


rm -rf build/ dist/ eeisp.egg-info/
python setup.py bdist_wheel
twine upload --repository testpypi dist/*

#twine upload --repository pypi dist/*
