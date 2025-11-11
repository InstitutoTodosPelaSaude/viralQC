import github
import os
import base64
import argparse
import requests


def download_dataset(repo_path: str, dataset_name: str, output_dir: str) -> None:
    """
    Download entire folder from GitHub with minimal API calls

    Args:
        repo_path (str): "owner/repo"
        dataset_name (str): Path to folder in repository
        output_dir (str): Local directory to save files
    """
    g = github.Github()
    repo = g.get_repo(repo_path)
    os.makedirs(output_dir, exist_ok=True)

    contents = repo.get_contents(dataset_name)
    for file_content in contents:
        local_path = os.path.join(output_dir, file_content.name)

        if file_content.encoding == "base64" and file_content.content:
            content = base64.b64decode(file_content.content)
        elif (
            file_content.encoding == "none"
            and hasattr(file_content, "_json_data")
            and "content" in file_content._json_data
        ):
            raw_content = file_content._json_data["content"]
            if raw_content:
                content = base64.b64decode(raw_content)
            else:
                raise Exception("No content in JSON data")
        else:
            file_data = repo.get_contents(file_content.path)
            if file_data.encoding == "base64" and file_data.content:
                content = base64.b64decode(file_data.content)
            else:
                raw_url = f"https://raw.githubusercontent.com/{repo_path}/main/{file_content.path}"
                response = requests.get(raw_url)
                if response.status_code == 200:
                    content = response.content

        with open(local_path, "wb") as f:
            f.write(content)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get nextclade dataset from a github repository."
    )
    parser.add_argument(
        "--repository",
        type=str,
        required=True,
        help="Repository name in the format 'owner/repo'.",
    )
    parser.add_argument(
        "--dataset-path",
        type=str,
        required=True,
        help="Path on repository to the dataset folder.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Path to store downloaded dataset.",
    )
    args = parser.parse_args()

    download_dataset(
        repo_path=args.repository,
        dataset_name=args.dataset_path,
        output_dir=args.output_dir,
    )
