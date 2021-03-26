
import pandas as pd

class SampleSheet():

    def __init__(self, path):
        sheet = pd.read_csv(path, sep='\t', header=0, index_col=0)
        assert sheet.index.is_unique
        self.sheet = sheet

    def getIds(self, filter=None):
        '''
        retrieve IDs (the index) of samples in the sample sheet.
        We can pass a tuple "filter" to retrieve only those for which a column (filter[0]) matches a certain value (filter[1]).

        returns a list of IDs.
        '''
        if filter:
            return self.sheet[self.sheet[filter[0]] == filter[1]].index.values.tolist()
        else:
            return self.sheet.index.values.tolist()

    def get(self, col, id=None):
        '''
        retrieve a column in the samplesheet, optionally only for specific samples (IDs, i.e. rows).
        '''
        if id:
            return self.sheet.loc[id, col]
        else:
            return self.sheet[col]

    def get_unique(self, col):
        '''
        get unique values of a column in the samplesheet.
        '''
        return self.sheet[col].unique().tolist()

    def __len__(self):
        return len(self.sheet.index)
    