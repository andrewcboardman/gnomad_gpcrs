init:
    pip install -r requirements.txt

test:
    python custom_constraint_analysis.py --test

standard:
	python custom_constraint_analysis.py --standard

.PHONY: init test standard