#!/usr/bin/env python

"""Utility classes to parse SAMPL submissions."""


# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import os
import io
import glob

import numpy as np
import pandas as pd

from .stats import mean_confidence_interval


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def load_submissions(submission_cls, directory_path, user_map):
    submissions = []
    for file_path in glob.glob(os.path.join(directory_path, '*.txt')):
        try:
            submission = submission_cls(file_path, user_map)
        except IgnoredSubmissionError:
            continue
        submissions.append(submission)
    return submissions


# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def plot_correlation(x, y, data, title=None, color=None, shaded_area_color=None, hue=None, ax=None):
    import seaborn as sns
    from matplotlib import pyplot as plt

    # Extract only free energies.
    values = data[[x, y]]

    # Find extreme values to make axes equal.
    # Generally plot between -20 and 0, and then
    # set it to the next number divisible by 5.
    min_limit = min(-20, np.floor(min(values.min())/5) * 5)
    max_limit = max(0, np.ceil(max(values.max())/5) * 5)
    axes_limits = np.array([min_limit, max_limit])

    if hue is None:
        grid = sns.jointplot(x=x, y=y, data=data,
                             kind='reg', joint_kws={'ci': None}, stat_func=None,
                             xlim=axes_limits, ylim=axes_limits, color=color)
        ax = grid.ax_joint
        grid.fig.subplots_adjust(top=0.95)
        grid.fig.suptitle(title)
    else:
        unique_hues = sorted(data[hue].unique())
        if ax is None:
            fig, ax = plt.subplots()
        # Set axes limits and ratio.
        ax.set_xlim(axes_limits)
        ax.set_ylim(axes_limits)
        ax.set_aspect('equal', 'box')
        # If a single color is passed, transform it into a palette.
        if not isinstance(color, list):
            color = [color for _ in range(len(unique_hues)+1)]
        # Add regression line single hue and all.
        for value, c in zip(unique_hues, color):
            sns.regplot(x=x, y=y, data=data[data[hue] == value], ci=0, label=value,
                        scatter=True, color=c, line_kws={'alpha': 0.5}, ax=ax)
        # Plot regression line for all the hues together.
        sns.regplot(x=x, y=y, data=data, ci=0, label='All',
                    scatter=False, color=color[len(unique_hues)],
                    line_kws={'alpha': 0.7}, ax=ax)
        ax.legend(loc='upper left')
        ax.set_title(title)

    # Add diagonal line.
    ax.plot(axes_limits, axes_limits, ls='--', c='black', alpha=0.8, lw=0.8)

    # Add shaded area for 1.5 kcal/mol error.
    if shaded_area_color is None:
        shaded_area_color = sns.color_palette('BuGn_r')[2]
    ax.fill_between(axes_limits, axes_limits - 1.5, axes_limits + 1.5, alpha=0.3,
                    color=shaded_area_color)


# =============================================================================
# MAIN SAMPL SUBMISSION CLASS
# =============================================================================

class IgnoredSubmissionError(Exception):
    """Exception used to signal a submission that must be ignored."""
    pass


class BadFormatError(Exception):
    """Exception used to signal a submission with unexpected formatting."""
    pass


class SamplSubmission:
    """A generic SAMPL submission.

    Parameters
    ----------
    file_path : str
        The path to the submission file.

    Raises
    ------
    IgnoredSubmission
        If the submission ID is among the ignored submissions.

    """
    # The D3R challenge IDs that are handled by this class.
    CHALLENGE_IDS = {}

    # The IDs of the submissions used for testing the validation.
    TEST_SUBMISSIONS = {}

    # Section of the submission file.
    SECTIONS = {}

    # Sections in CSV format with kwargs to pass to pandas.read_csv().
    CSV_SECTIONS = {}

    def __init__(self, file_path, user_map):
        file_name = os.path.splitext(os.path.basename(file_path))[0]
        file_data = file_name.split('-')

        # Check if this is a deleted submission.
        if file_data[0] == 'DELETED':
            raise IgnoredSubmissionError('This submission was deleted.')

        # Check if this is a test submission.
        self.receipt_id = file_data[0]
        if self.receipt_id in self.TEST_SUBMISSIONS:
            raise IgnoredSubmissionError('This submission has been used for tests.')

        # Check this is the correct challenge.
        self.challenge_id = int(file_data[1])
        assert self.challenge_id in self.CHALLENGE_IDS

        # Store user map information.
        if user_map is not None:
            user_map_record = user_map[user_map.receipt_id == self.receipt_id]
            assert len(user_map_record) == 1
            user_map_record = user_map_record.iloc[0]

            self.id = user_map_record.id
            self.participant = user_map_record.firstname + ' ' + user_map_record.lastname
            self.participant_id = user_map_record.uid
            self.participant_email = user_map_record.email
            assert self.challenge_id == user_map_record.component
        else:
            self.id = None
            self.participant = None
            self.participant_id = None
            self.participant_email = None

    @classmethod
    def _read_lines(cls, file_path):
        """Generator to read the file and discard blank lines and comments."""
        with open(file_path, 'r', encoding='utf-8-sig') as f:
            for line in f:
                # Strip whitespaces.
                line = line.strip()
                # Don't return blank lines and comments.
                if line != '' and line[0] != '#':
                    yield line

    @classmethod
    def _load_sections(cls, file_path):
        """Load the data in the file and separate it by sections."""
        sections = {}
        current_section = None
        for line in cls._read_lines(file_path):
            # Check if this is a new section.
            if line[:-1] in cls.SECTIONS:
                current_section = line[:-1]
            else:
                try:
                    sections[current_section].append(line)
                except KeyError:
                    sections[current_section] = [line]

        # Check that all the sections have been loaded.
        found_sections = set(sections.keys())
        if found_sections != cls.SECTIONS:
            raise BadFormatError('Missing sections: {}.'.format(cls.SECTIONS - found_sections))

        # Create a Pandas dataframe from the CSV format.
        for section_name, pandas_kwargs in cls.CSV_SECTIONS.items():
            csv_str = io.StringIO('\n'.join(sections[section_name]))
            section = pd.read_csv(csv_str, skipinitialspace=True, **pandas_kwargs)
            sections[section_name] = section
        return sections

    @classmethod
    def _create_comparison_dataframe(cls, column_name, submission_data, experimental_data):
        """Create a single dataframe with submission and experimental data."""
        # Filter only the systems IDs in this submissions.
        experimental_data = experimental_data[experimental_data.index.isin(submission_data.index)]
        # Fix the names of the columns for labelling.
        submission_series = submission_data[column_name]
        submission_series.name += ' (calc)'
        experimental_series = experimental_data[column_name]
        experimental_series.name += ' (expt)'
        # Concatenate the two columns into a single dataframe.
        return pd.concat([experimental_series, submission_series], axis=1)

