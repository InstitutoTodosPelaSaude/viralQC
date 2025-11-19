import os
import argparse
import zipfile
import requests
import tempfile


def download_dataset(repo_path: str, dataset_name: str, output_dir: str) -> None:
    """
    Download a GitHub folder without using the GitHub API.

    Args:
        repo_path (str): "owner/repo"
        dataset_name (str): Path to folder in repository
        output_dir (str): Local directory to save files
    """
    owner = repo_path.split("/")[0]
    repo = repo_path.split("/")[1]

    zip_url = f"https://codeload.github.com/{owner}/{repo}/zip/refs/heads/main"
    resp = requests.get(zip_url, stream=True)

    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        for chunk in resp.iter_content(1024 * 64):
            tmp.write(chunk)
        zip_path = tmp.name

    zip_prefix = f"{repo}-main/{dataset_name}/"
    os.makedirs(output_dir, exist_ok=True)

    with zipfile.ZipFile(zip_path) as zf:
        for name in zf.namelist():
            if (
                name.startswith(zip_prefix)
                and name != zip_prefix
                and not name.endswith("/")
            ):
                relative_name = name[len(zip_prefix) :]
                # ignore files in subfolders (non-recursive)
                if "/" in relative_name:
                    continue

                out_file = os.path.join(output_dir, relative_name)
                with zf.open(name) as src, open(out_file, "wb") as dst:
                    dst.write(src.read())

    os.remove(zip_path)


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
