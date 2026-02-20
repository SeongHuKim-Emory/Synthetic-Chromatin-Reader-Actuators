# git_push.ps1

# Ensure PowerShell is running the script as an executable
# Set permissions with: Set-ExecutionPolicy RemoteSigned -Scope CurrentUser

# Retrieve and display the remote repository address
$remoteRepo = git remote get-url origin
if (-not $remoteRepo) {
    Write-Host "No remote repository configured. Exiting."
    exit 1
}
Write-Host "Remote repository: $remoteRepo"

# Prompt the user for input
$COMMIT_MESSAGE = Read-Host "Enter the Commit Message"
$TAG = Read-Host "Enter the Tag"

# Check if inputs are empty
if (-not $COMMIT_MESSAGE -or -not $TAG) {
    Write-Host "Both Commit Message and Tag are required. Exiting."
    exit 1
}

# Execute Git commands
git pull origin main
git add .
git commit -m $COMMIT_MESSAGE
git branch -M main
git push -u origin main
git tag $TAG
git push origin $TAG

# Completion message
Write-Host "Git operations completed successfully."

