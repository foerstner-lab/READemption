test:
	pytest

coverage:
	coverage run tests/test_all.py
	coverage report

package:
	pandoc -o README.rst README.md
	python3 setup.py sdist bdist_wheel
	rm -rf READemption.egg-info
	ls dist/*
	rm README.rst

package_to_pypi:
	python3 -m twine upload  --verbose  dist/*

html_doc:
	cd docs && make html && cd ..

show_html_docs:
	firefox docs/build/html/index.html &

pylint:
	pylint bin/reademption reademptionlib/* tests/*

new_release:
	@echo "* Create/checkout a release branch"
	@echo "  git branch release_v0.4.X"
	@echo "  git checkout release_v0.4.X"
	@echo "* Change version and date in bin/reademption"
	@echo "* Change version in setup.py"
	@echo "* Change version, release and date in docs/source/conf.py"
	@echo "* Add description to CHANGELOG.txt"
	@echo "* Create new docs"
	@echo "  make html_doc" # was this only needed for pythonhosted?
	@echo "* Create package"
	@echo " make package"
	@echo "*delete old package"
	@echo "rm dist/READemption-X.0.0-py3-none-any.whl"  
	@echo "rm dist/READemption-X.0.0.tar.gz"
	@echo "* Upload package to pypi"
	@echo " make package_to_pypi"
	@echo "* Test doc creation"
	@echo "* git add CHANGELOG.txt bin/reademption docs/source/conf.py setup.py"
	@echo " Also add the changed html pages of the docu"
	@echo "* Commit changes e.g. 'git commit -m \"Set version to 0.4.X\"'"
	@echo "* Tag the commit e.g. 'git tag -a v0.4.X -m \"version v0.4.X\"'"
	@echo "* Merge release into dev and master"
	@echo "* Push it to github: git push"
	@echo "* Generate a new release based on this tag at"
	@echo "  https://github.com/foerstner-lab/READemption/releases/new"
	@echo "Update the documentation on read the docs"
	      " -  login and hit 'build'"
	# @echo "* Upload new docs using 'make upload_doc'"

# Todo: docker_images: docker build -t konradfoerstner/reademption:0.4.1 -t konradfoerstner/reademption:latest .
