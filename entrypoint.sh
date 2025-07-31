#!/bin/bash
set -e

# Ensure the use of global pip and python, without virtual environment
python manage.py migrate
python manage.py runserver 0.0.0.0:8000
