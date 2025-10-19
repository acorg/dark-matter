import nox

nox.options.default_venv_backend = "uv"
nox.options.reuse_existing_virtualenvs = True


@nox.session(python=["3.10", "3.11", "3.12"])
def tests(session):
    session.install(".", "pytest")
    session.run("pytest")


if __name__ == "__main__":
     nox.main()
