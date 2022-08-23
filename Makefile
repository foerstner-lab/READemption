test:
	pytest

coverage:
	coverage run tests/test_all.py
	coverage report

coverage_html:
	# pytest-cov needs to be installed:
	# $ pip install pytest-cov
	pytest --cov=reademptionlib --cov-report html --cov-context=test

coverage_batch:
	# install coverage-badge
	# pip install coverage-badge
	# run coverage
	coverage run tests/test_all.py
	# create the badge
	coverage-badge -o coverage.svg
	# add and commit the badge to GitHub

package:
	rm -r build
	pandoc -o README.rst README.md
	python3 setup.py sdist bdist_wheel
	rm -rf READemption.egg-info
	ls dist/*
	rm README.rst

package_to_pypi:
	python3 -m twine upload  --verbose  dist/*

html_doc:
	# install sphinx argparse extension
	# pip install sphinx-argparse
	# install rtd theme
	# pip install sphinx-rtd-theme
	cd docs && make html && cd ..

show_html_docs:
	firefox docs/build/html/index.html &

pylint:
	pylint bin/reademption reademptionlib/* tests/*

conda_package:
	# create conda meta.yaml from pypi package
	conda skeleton pypi reademption

new_release:
	@echo "* Create/checkout a release branch"
	@echo "  git checkout dev"
	@echo "* Change version and date in bin/reademption"
	@echo "* Change version in setup.py"
	@echo "* Change requirements (install_requires) in setup.py"
	@echo "* Change version, release and date in docs/source/conf.py"
	@echo "* Add description to CHANGELOG.txt"
	@echo "* Refresh the two reademption links in docs and docs/source."
	      "*  this is need for building the docu on 'read the docs'. See last point of this list."
	@echo "rm docs/reademption && rm docs/source/reademption && cd docs && ln -s ../bin/reademption reademption && cd source && ln -s ../../bin/reademption reademption && cd ../.."
	@echo "* Create new docs"
	@echo "* if making the html documentation fails, try to install the necessary programs:"
		  "* $ pip install sphinx sphinx_rtd_theme"
		  "* $ pip install sphinx-argparse"
	@echo "  make html_doc"
	@echo "* Create package"
	@echo " make package"
	@echo "*delete old package"
	@echo "rm dist/READemption-X.0.0-py3-none-any.whl"  
	@echo "rm dist/READemption-X.0.0.tar.gz"
	@echo "* Upload package to pypi"
	@echo "* If necessary install twine e.g.: "
	@echo "* conda install -c conda-forge twine "
	@echo " make package_to_pypi"
	@echo "* Test doc creation"
	@echo "* Update version and date in CITATION.cff"
	@echo "* git add CHANGELOG.txt bin/reademption docs/source/conf.py setup.py CITATION.cff"
	@echo " Also add the changed html pages of the docu e.g.: git add -u    (adds all modified files to repo)"
	@echo "* Commit changes e.g. 'git commit -m \"Set version to 0.4.X\"'"
	@echo "* Merge release into dev and main"
	@echo " git checkout main"
	@echo "* delete docs/build to avoid conflicts with new binary doc files"
	@echo "git rm -r docs/build"
	@echo "git commit -m 'Deleted old docs build'"
	@echo " git merge dev"
	@echo "Solve conflicts from merging (CHANGELOG.txt bin/reademption docs/source/conf.py setup.py CITATION.cff)"
	@echo "* Tag the commit e.g. 'git tag -a v0.4.X -m \"version v0.4.X\"'"
	@echo "* Push it to github: git push"
	@echo "* Generate a new release based on this tag at"
	@echo "  https://github.com/foerstner-lab/READemption/releases/new"
	@echo "The documenation should update automatically after pushing to github via a .readthedocs.yml"
	@echo "If the automatic update did not work go to read the docs"
	      " -  login and hit 'build'"
	@echo "*create conda package"
	@echo "Change version in conda_package/reademption/meta.yaml"
	@echo "conda-build --python 3.9 conda_package/reademption"
	@echo "if necessary install anaconda client"
	@echo " conda install anaconda-client"
	@echo "* login to anaconda client"
	@echo " anaconda login"
	@echo "*upload anaconda package"
	@echo "anaconda upload /home/till/anaconda3/conda-bld/linux-64/reademption-2.0.0-py39_0.tar.bz2"
	@echo "*Create Docker image with tags latest and the current version:"
	@echo " docker build -f Dockerfile -t tillsauerwein/reademption:2.0.0 -t tillsauerwein/reademption:latest ."
	@echo "*push all docker images to dockerhub:"
	@echo " docker push -a tillsauerwein/reademption"
	# @echo "* Upload new docs using 'make upload_doc'"

