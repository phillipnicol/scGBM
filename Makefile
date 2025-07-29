.PHONY: pull push

pull:
	git pull
	Rscript -e "renv::restore()"
	R CMD INSTALL .

push:
	Rscript -e "renv::snapshot()"
	Rscript -e "roxygen2::roxygenize()"
	git add .
	git commit -m "Auto-commit from makefile"
	git push
