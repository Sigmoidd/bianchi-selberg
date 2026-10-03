"""Compatibility command for the d=2 certificate; uses the uniform CLI."""
if __package__:
    from .group_certificate import main
else:
    from group_certificate import main

if __name__ == "__main__":
    main(default_field=2)
