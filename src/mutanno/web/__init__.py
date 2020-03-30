#!/usr/bin/env python
"""Django's command-line utility for administrative tasks."""
import os
import sys
# import django


def runserver(opt):
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'mutanno.web.settings')
    # django.setup()

    try:
        from django.core.management import execute_from_command_line
    except ImportError as exc:
        raise ImportError(
            "Couldn't import Django. Are you sure it's installed and "
            "available on your PYTHONPATH environment variable? Did you "
            "forget to activate a virtual environment?"
        ) from exc
    argvlist = ["", "runserver", opt['port']]
    execute_from_command_line(argvlist)
