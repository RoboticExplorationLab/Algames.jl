# Adding Documentation
This section details how to set up a documentation page using [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl). This is a more concise version of the instructions included
in the documentation for that package.

1. In your repo, create a `docs/` directory, and a `docs/src/` directory that will conntain
    all the source `.md` files for your documentation.

2. Create a `docs/make.jl` file. This file is responsible for building and deploying your
    documentation. Here is the basic starting point:

        using Documenter
        using Algames  # your package name here

        makedocs(
            sitename = "Algames",  # your package name here
            format = Documenter.HTML(prettyurls = false),  # optional
            pages = [
                "Introduction" => "index.md"
            ]
        )

        # Documenter can also automatically deploy documentation to gh-pages.
        # See "Hosting Documentation" and deploydocs() in the Documenter manual
        # for more information.
        deploydocs(
            repo = "github.com/bjack205/Algames.jl.git",
        )

3. Add documentation files to `docs/src`. Once the files are in `docs/src`, add them to
    the `makedocs` command.

4. Add Documentation dependencies. Nearly identical to the tests, we need to add any
    dependencies we use to build the documentation, which obviously must include Documenter.jl.
    Activate the `docs/` directory and add Documenter
    ```
    julia> ] activate docs
    (docs) pkg> add Documenter
    ```
    Then add a `[compat]` entry for Documenter.

4. Add deploy keys for your repo. Install [DocumenterTools.jl](https://github.com/JuliaDocs/DocumenterTools.jl) and enter the following into your REPL
    ```
    using DocumenterTools
    using Algames                     # your package name here
    DocumenterTools.genkeys(Algames)  # your package name here
    ```
    Copy the first public key (starts with `ssh-rsa` and ends with ` Documenter`).
    Go to your repository settings in GitHub and select "Deploy Keys". Add the deploy key,
    using `documenter` as the name.
    ![deploy_key](images/deploy_key.png)

    Then copy the very long environment variable and save it as the `DOCUMENTER_KEY` secret
    on GitHub:
    ![doc_key](images/doc_key.png)

5. Add a GitHub Action to build your documentation. Create a new GitHub action
    (called `Documenter.yml`) and paste the code found
    [here](https://github.com/bjack205/Algames.jl/blob/master/.github/workflows/Documenter.yml):

6. Add Documentation badge to README. Add the following line to the top of the file,
    replacing the user/organize and repo names in the url:
    `[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://bjack205.github.io/Algames.jl/dev)`

## Writing Docstrings
As stated in the [Julia manual](https://docs.julialang.org/en/v1/manual/documentation/),
start the docstring with the signature, which is indented with 4 spaces so it prints as
Julia code. You can use normal markdown syntax, such as headings, to make your docstrings
look since and stay organized. See the above example and the
[Documenter.jl docs](https://juliadocs.github.io/Documenter.jl/stable/man/latex/) on how
to include ``\LaTeX`` math into the docstrings.



## Building documentation locally
You can build the documentation locally by only running the `makedocs` function, and
disabling `prettyurls`. It's common to update docstrings in your code and want these
changes reflected in your build. After making a change to the docstring, you need to
"rebuild" the docstrings by executing the whole file, easily done with `CTRL-SHIFT-RETURN`
in Juno. You can then rebuild the the docs (using `CTRL-RETURN` in Juno) and the docstrings
will be updated.
