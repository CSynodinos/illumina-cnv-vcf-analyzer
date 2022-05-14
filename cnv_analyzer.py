#!/usr/bin/env python3
from __future__ import annotations

import sys
from inspect import getfullargspec
from re import search, sub, compile, split
import pandas as pd


class NoCnvFoundError(Exception):
    """Custom exception class raised when no CNV's are present in a .vcf file."""

    __module__ = 'builtins'

    def __init__(self, *args) -> None:
        if args:
            self.errmessage = args[0]
        else:
            self.errmessage = None

    def __repr__(self) -> str:
        if self.errmessage:
            return '{0} '.format(self.errmessage)
        else:
            return 'NoCnvError has been raised'

class vcf_parser:

    cnvs = ('<DUP>', '<DEL>')

    def __init__(self, vcf_fl: str, disp_info = None, find_only_dups = False, find_only_dels = False,
                out = 'output.txt') -> None:
        self.vcf_fl = vcf_fl
        self.disp_info = disp_info
        self.find_only_dups = find_only_dups
        self.find_only_dels = find_only_dels
        self.out = out

    @classmethod
    def __repr__(cls) -> str:
        params = getfullargspec(__class__).args
        params.remove("self")
        return params

    @classmethod
    def __dir__(cls, only_added = False) -> list:
        """Display function attributes.
        Args:
            * `only_added` (bool, optional): Choose whether to display only the specified attributes. Defaults to False.
        Returns:
            list: List of attributes.
        """

        all_att = list(cls.__dict__.keys())
        if not only_added:
            return all_att
        else:
            default_atts = ['__module__', '__doc__', '__dict__', '__weakref__']
            all_att = [x for x in all_att if x not in default_atts]
            return all_att

    @classmethod
    def raise_cnv_err(cls, f):
        """Raises NoCnvFoundError when called.

        Args:
            * `f` (str): Name of file used in error message.

        Raises:
            `NoCnvFoundError`: No cnv's are present in the .vcf file.
        """
        raise NoCnvFoundError(f"Unable to find any cnv's in file: {f}. File must contain at least one of the following cnv's: {', '.join(cls.cnvs)}")

    def cnv_line_finder(self) -> tuple[list, str]:
        """Finds all the lines in a .vcf file that contain cnv's.

        Raises:
            * `RuntimeError`: Raised when the script is unable to find the first chromosome entry. Potential compatibility issue with the .vcf file.
            * `NoCnvFoundError`: No cnv's are present in the .vcf file.

        Returns:
            tuple[list, str]: List of all the lines with cnv entries, string equal to the name of the sample.
        """

        with open(self.vcf_fl, "r") as vcf:
            out_lst = []
            counter = 0
            for line in vcf:
                line = line.rstrip()
                if "#CHROM" in line:
                    sample = line.split("\t")[-1]
                if not line[0:3] == "chr":
                    continue
                counter += 1
                if "<DUP>" in line or "<DEL>" in line:
                    out_lst.append(line)

        if counter == 0:
            raise RuntimeError("Unable to locate the first chromosome entry. Every chromosome entry must start with the chr character.")
        if not len(out_lst) > 0:
            return self.raise_cnv_err(f = self.vcf_fl)

        return out_lst, sample

    @staticmethod
    def _add_to_df(iter: str, d: pd.DataFrame) -> bool:
        row_no_tabs = iter.replace('\t', " ").split(" ")
        d.loc[len(d)] = row_no_tabs
        return True

    def _info_parser(self):
        lst_of_cnv_lines, sample = self.cnv_line_finder()
        indv_entities = {}  # fill with whatever the user asks for. Only use for one call that depends on specific start. give back the full line. Each element is 1 cell. 
        df = pd.DataFrame([[0,0,0,0,0,0,0,0,0,0]], columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample])
        for i in lst_of_cnv_lines:
            chrom = compile(r"^[^\t]*")
            chrom = search(chrom, i).group(0)
            startpos = split("\t+", i)
            start = []
            counter = 0
            for pos in startpos:
                if pos.isdigit() and counter == 0:  # Counter is only equal to 0 when the start has not been found, once it gets found this block doesn't run again.
                    start.append(pos)
                    start = "".join(start)
                    counter += 1

                if "DRAGEN:" in pos:
                    start_end = sub("[^0-9]", ", ", pos)
                    start_end = start_end.replace(",","").strip()
                    whitespace_sep = compile(r"^[^ ]*")
                    start_end = start_end.split(" ", 1)[1]
                    start_ref = search(whitespace_sep, start_end).group(0)  # start and end will be shown at the end for specific calls.
                    end = sub(start_ref, "", start_end).strip()
                    if not self.disp_info == None:  # User did not request specific info:
                        if self.disp_info.isdigit() and (self.disp_info == start or self.disp_info == end):
                            indv_entities["sample"] = sample
                            indv_entities["chromosome"] = chrom
                            indv_entities["start"] = start
                            indv_entities["end"] = end
                            if "<DUP>" in i:
                                indv_entities["cnv_type"] = "<DUP>"
                            else:
                                indv_entities["cnv_type"] = "<DEL>"
                            qc_score = i.split(end, 1)[1]
                            qc_score = qc_score.split("SVLEN", 1)[0]
                            qc_score = sub("[^0-9]", "", qc_score).strip()
                            indv_entities["score"] = qc_score
                            if "cnvQual" in i:
                                indv_entities["filter"] = "cnvQual"
                            elif "PASS" in i:
                                indv_entities["filter"] = "PASS"

                        if not self.disp_info.isdigit():
                            raise TypeError(f"-inf argument must be an integer. {self.disp_info} is not an integer.")
                else:
                    continue

            if self.find_only_dups:
                if "<DUP>" in i:
                    self._add_to_df(iter = i, d = df)
            elif self.find_only_dels:
                if "<DEL>" in i:
                    self._add_to_df(iter = i, d = df)
            elif not (self.find_only_dups and self.find_only_dels):
                if "<DUP>" in i or "<DEL>" in i:
                    self._add_to_df(iter = i, d = df)

        df = df.iloc[1: , :]    # Removes the first place-holder row.
        if len(indv_entities) > 0:
            return df, indv_entities
        else:
            return df

    def _stats_writer(self):
        if self.disp_info:
            df, single_entry = self._info_parser()
        else:
            single_entry = None
            df = self._info_parser()

        # To be used in description.
        num_of_df_entries = len(df.index)   # Number of CNV's

        ## cnv's
        all_cnvs = df["ALT"].value_counts()
        if self.find_only_dups:
            ent_index = ["DUPS"]
            all_cnvs.index = ent_index
            num_of_df_dups = str(all_cnvs.get(key = "DUPS"))
            num_of_df_dels = None
            qual_vals = df["QUAL"].tolist()
            qual_vals = [int(x) for x in qual_vals]
            max_dups = max(qual_vals)
            max_dels = None
            min_dups = min(qual_vals)
            min_dels = None
        elif self.find_only_dels:
            ent_index = ["DELS"]
            all_cnvs.index = ent_index
            num_of_df_dups = None
            num_of_df_dels = str(all_cnvs.get(key = "DELS"))
            qual_vals = df["QUAL"].tolist()
            qual_vals = [int(x) for x in qual_vals]
            max_dups = None
            max_dels = max(qual_vals)
            min_dups = None
            min_dels = min(qual_vals)
        else:
            ent_index = ["DUPS", "DELS"]
            all_cnvs.index = ent_index
            num_of_df_dups = str(all_cnvs.get(key = "DUPS"))
            num_of_df_dels = str(all_cnvs.get(key = "DELS"))
            print(df['ALT',])
        ## Scores
        if not "DUPS" in set(all_cnvs):
            dups_score = None
            dels_score = df["QUAL"].value_counts().describe()
        if not "DELS" in set(all_cnvs):
            dups_score = df["QUAL"].value_counts().describe()
            dels_score = None
        if "<DEL>" and "<DUP>" in set(df["ALT"]):
            df_scores = df.loc[:, df.columns.intersection(['ALT','QUAL'])]
            df_dups = df_scores.loc[df['ALT'] == "<DUP>"].drop(columns = 'ALT')
            df_dels = df_scores.loc[df['ALT'] == "<DEL>"].drop(columns = 'ALT')
            dups_score = df_dups["QUAL"].value_counts().describe()
            dels_score = df_dels["QUAL"].value_counts().describe()
        
        def __count_cnv(score: pd.Series):
            """Internal function to count the number of a specified cnv.

            Args:
                score (pd.Series): A column from a pandas Series.

            Returns:
                str: The number of occurances of a cnv rounded to 0 decimal places.
            """

            _count = list(score.loc[['count']])
            _count = ''.join([str(c) for c in _count])
            return str(int(float(_count)))

        if not dups_score.empty:
            unique_dups_count = __count_cnv(score = dups_score)
        if not dels_score.empty:
            unique_dels_count = __count_cnv(score = dels_score)

        with open(self.out, "w") as txt:
            txt.write('Copy Number Variant Analyzer')
            txt.write('\n')
            txt.write('This file was produced with the cnv_analyzer python tool.')
            txt.write('\n')
            txt.write('\n')
            txt.write('Summary of analysis:')
            txt.write('\n')
            txt.write('\n')

            txt.write("The analyzer looked for the following types of cnv's: ")
            if self.find_only_dups:
                txt.write('DUP')
            if self.find_only_dels:
                txt.write('DEL')
            if not (self.find_only_dups and self.find_only_dels):
                txt.write('DUP, DEL')

            txt.write('\n')
            txt.write(f"Total number of CNV's: {num_of_df_entries}")
            txt.write('\n')
            if num_of_df_dups == None:
                txt.write(f'Number of deletion entries found: {num_of_df_dels}')
                txt.write('\n')
                txt.write(f'Number of unique deletion entries: {unique_dels_count}')
            if num_of_df_dels == None:
                txt.write(f'Number of duplication entries: {num_of_df_dups}')
                txt.write('\n')
                txt.write(f'Number of unique deletion entries: {unique_dups_count}')
            if not (num_of_df_dups and num_of_df_dels) == None:
                txt.write(f'Number of duplication entries: {num_of_df_dups}')
                txt.write('\n')
                txt.write(f'Number of deletion entries: {num_of_df_dels}')
                txt.write('\n')
                txt.write(f'Number of unique duplication entries: {unique_dups_count}')
                txt.write('\n')
                txt.write(f'Number of unique deletion entries: {unique_dels_count}')


def main():
    vcf_parser(vcf_fl="test_files/a_sample.cnv.vcf")._stats_writer()

if __name__ == "__main__":
    main()