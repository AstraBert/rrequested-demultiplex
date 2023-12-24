import os


class FileSystem:
    def __init__(self):
        pass

    def makedir_orchange(self, path):
        """
        Creates a directory at the specified path if it doesn't exist.
        Returns the provided path.

        Parameters:
        - path (str): The directory path to be created.

        Returns:
        - str: The provided path.
        """
        try:
            os.mkdir(path)
            return path
        except FileExistsError:
            return path

    def get_base_dir(self, path):
        """
        Extracts information about the base directory, base name (without extension),
        and the base directory of the base directory from the given file path.

        Parameters:
        - path (str): The file path for which information is to be extracted.

        Returns:
        - tuple: A tuple containing the base directory, base name, and the base directory of the base directory.
        """
        base, ext = os.path.splitext(path)
        asedir = base.split("/")
        b = "/"
        basedir = b.join(asedir[:len(asedir) - 1])
        base_basedir = b.join(asedir[:len(asedir) - 2])
        basename = asedir[len(asedir) - 1]
        return basedir, basename, base_basedir

