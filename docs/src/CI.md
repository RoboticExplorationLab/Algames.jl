# Setting Up CI
With a functional test suite, we are now ready to set up CI. Here we will demonstrate
how to do this with GitHub actions.

1. Open GitHub Actions panel in GitHub:
    ![GHA](images/GHA_button.png)

2. Click "set up a workflow yourself"
    ![GHA_setup](images/GHA_setup.png)

3. Rename the file to `CI.yml`
    ![GHA_rename](images/GHA_rename.png)

4. Paste this [code](https://github.com/bjack205/Algames.jl/blob/master/.github/workflows/CI.yml).
    Which runs Julia v1.3 (under `jobs/test/strategy/matrix/version`) on Ubuntu, Mac, and
    Windows (under `jobs/test/strategy/matrix/os`). The steps setup Julia, run the tests,
    and then upload the code coverage results to codecov. This will run on any pull request
    and on any push to the master branch.

5. Commit the file and GitHub will automatically start running the tests

6. Add the badge: Go back to the GitHub Actions pane and select the currently running
   CI action and select "Create Status Badge", copy the badge link and paste it into the
   top of your README file.
   ![GHA_badge](images/GHA_badge.png)
