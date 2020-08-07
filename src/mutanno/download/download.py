import urllib.request
import os
from tqdm import tqdm

from ..util import file_util
from .. import preprocess


class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_file(url, output_path):
    with DownloadProgressBar(unit='B', unit_scale=True,
                             miniters=1, desc=url.split('/')[-1]) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)


def get_download_filename(dirpath, refversion, source_name, file_version, url):
    fname = url.split('/')[-1]
    outfile = os.path.join(dirpath, 'downloaded', source_name, refversion, file_version, fname)
    return outfile


class Downloader():
    def __init__(self, opt):
        self.opt = opt
        self.sourcelist = {}

    def download_sourcelist_file(self):
        pass

    def load_sourcelist_file(self):
        sourcelist_filepath = file_util.getDataPath('sourcelist.json')
        self.sourcelist = file_util.load_json(sourcelist_filepath)

    def download_source_file(self, source_name, file_version, url):
        outfile = get_download_filename(self.opt.dir, self.opt.refversion, source_name, file_version, url)
        flag_download = True
        if file_util.is_exist(outfile):
            msg = "The file (" + outfile + ") exists.\n Do you want to overwrite it. (y/n)?"
            val = input(msg)
            if val.upper() != "Y":
                flag_download = False
        if flag_download:
            file_util.check_dir(outfile)
            download_file(url, outfile)
        return outfile

    def preprocess_source_file(self, function_name, rawfiles, source_name, file_version):
        run_function = getattr(preprocess, function_name)
        data = {}
        data['source_name'] = source_name
        data['refversion'] = self.opt.refversion
        data['rawfiles'] = rawfiles
        data['file_version'] = file_version
        data['dirpath'] = self.opt.dir
        data['outfile_title'] = '_'.join([source_name, self.opt.refversion, file_version])
        result = run_function(data)
        return result

    def download_and_preprocess_sourcefile(self):
        for s1 in self.sourcelist['sources'][self.opt.refversion]:
            if self.opt.source == "all" or self.opt.source == s1['name']:
                if self.opt.version == "latest":
                    selected_version = s1['latest']
                elif self.opt.version != "":
                    selected_version = self.opt.version

                for v1 in s1['versions']:
                    if v1['version'] == selected_version:
                        outfiles = []
                        for url in v1['urls']:
                            outfiles.append(self.download_source_file(s1['name'], v1['version'], url))
                        self.preprocess_source_file(s1['preprocess_function'], outfiles, s1['name'], v1['version'])

    def run(self):
        self.download_sourcelist_file()
        self.load_sourcelist_file()
        self.download_and_preprocess_sourcefile()
