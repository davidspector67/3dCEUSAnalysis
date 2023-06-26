unix: main.py
	source venv/bin/activate
	python main.py
	deactivate

windows: main.py
	call venv\scripts\activate.bat
	python main.py
	deactivate

build_unix: pyPackages.txt
	pip install virtualenv
	python -m venv venv
	source venv/bin/activate
	pip install -r pyPackages.txt
	deactivate

build_windows: pyPackages.txt
	pip install virtualenv
	python -m venv venv
	call \venv\scripts\activate.bat
	pip install -r pyPackages.txt
	deactivate