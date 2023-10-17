'''
Dendrogram annotations

This file holds the Annotation class, which builds the blueprint for plotting
dendrograms in iTOL.

Class
-----
Annotation
'''
from typing import Union, Tuple

class Annotation():
    '''
    Class that represents an annotation iTOL file.

    Attributes
    ----------
    separator : str
        Separator used to delimit elements of the data throughout the 
        annotation file. The choices are 'SPACE', 'TAB' and 'COMMA'. The 
        default is 'SPACE'.
    
    sep : str
        The character the separator represents. The choices are ' ', '\t' and 
        ','. The default is ' '.
    
    color : str
        Dataset color. The default is '#000' (black).
    
    margin : int
        Left margin, used to increase/decrease the spacing to the next dataset. 
        Can be negative, causing datasets to overlap. In plain words, if one 
        desires to add a colorstrip as well and binary labels on their tree, 
        this number controls the distance between them. The default is 0.
    
    show_internal : bool
        Show internal values. If set, values associated to internal nodes will 
        be displayed even if these nodes are not collapsed. It could cause 
        overlapping in the dataset display.The default is 0.
    
    color_branches : bool
        Choose whether one wishes to color the branches with the same color of
        the leaf nodes. The default is 0.
    
    strip_width : int
        Strip width in pixels. The default is 70.

    borther_width : int
        If set above 0, a black border of specified width (in pixels) will be 
        drawn around the color strip. The default is 10.

    border_color : str
        Color of the border lines around the color strip. The default is #00 
        (black).
    
    field_labels : Union[str, Tuple[str]]
        Name of each field. The default is ('C1', 'C4', 'substrate').

    field_colors : Union[str, Tuple[str]]
        Colot to fill the field shape. The default is ('#000', '#000', '#66f')

    field_shapes : Union[str, Tuple[str]]
        Shapes for each field column. Choices are {1: rectangle, 2: circle, 
        3: star, 4: right pointing triangle, 5: left pointing triangle}. The 
        default is (5, 4, 3).

    height_factor : int
        Symbol height factor. Default symbol height (not the one here) will be 
        slightly less than the available space between leaves, but you can set 
        a multiplication factor here to increase/decrease it (values from 0 to 
        1 will decrease it, values above 1 will increase it). The default is 
        85.
    
    symbol_spacing : int
        Increase/decrease the spacing between individual levels, when there is 
        more than one binary level defined. The default is 10.    

    Methods
    -------
    __init__
    colorstrip
    binary
    '''
    def __init__(
            # General attributes
            self,
            separator : str = 'SPACE',
            color : str = '#000',
            margin : int = 0,
            show_internal : bool= 0,

            # Colorstrip-specific attributes
            color_branches : bool = 0,
            strip_width : int = 70,
            border_width : int = 10,
            border_color : str = '#000',

            # Binary-specific attributes
            field_labels : Union[str, Tuple[str]] = ('C1', 'C4', 'substrate'),
            field_colors : Union[str, Tuple[str]] = ('#000', '#000', '#66f'),
            field_shapes : Union[str, Tuple[str]] = (4, 5, 3),
            height_factor : int = 85,
            symbol_spacing : int = 10

            ):
        
        separator_options = {'SPACE' : ' ', 'COMMA' : ',', 'TAB':'\t'}
        self.separator = separator
        self.sep = separator_options[separator]
        self.color = color
        self.margin = margin
        self.show_internal = show_internal
        self.color_branches = color_branches
        self.strip_width = strip_width
        self.border_width = border_width
        self.border_color = border_color
        self.field_labels = field_labels
        self.field_colors = field_colors
        self.field_shapes = field_shapes
        self.height_factor = height_factor
        self.symbol_spacing = symbol_spacing

    def colorstrip(self, data : dict[str, str]) -> str:
        '''
        Builds the colorstrip annotation file for an iTOL tree from the leaf
        names and their color.

        Parameters
        ----------
        data : dict[str, str]
            Leaf names as keys and color as values. In our specific situation, 
            the leaf names are usually UniProt IDs (e.g. A0A4P7DN87) and the
            color is one arbitrarily associated to the AA family the protein 
            belongs to (e.g. #1f77b4)

        Returns
        -------
        content : str
            Colorstrip annotation file content
        '''

        content = 'DATASET_COLORSTRIP\n'
        content += f'SEPARATOR {self.separator}\n'
        content += f'DATASET_LABEL{self.sep}colorstrip\n'
        content += f'COLOR{self.sep}{self.color}\n'
        content += f'COLOR_BRANCHES{self.sep}{self.color_branches}\n'
        content += f'STRIP_WIDTH{self.sep}{self.strip_width}\n'
        content += f'MARGIN {self.sep}{self.margin + 85}\n'
        content += f'BORDER_WIDTH{self.sep}{self.border_width}\n'
        content += f'BORDER_COLOR{self.sep}{self.border_color}\n'
        content += f'SHOW_INTERNAL{self.sep}{self.show_internal}\n'
        content += 'DATA\n'
        for branch in data:
            content += f'{branch}{self.sep}{data[branch]}\n'
        
        return content

    def binary(self, data : dict[str, str]) -> str:
        '''
        Builds the binary annotation file for an iTOL tree from the leaf
        names and their binary data.

        Parameters
        ----------
        data : dict[str, str]
            Leaf names as keys and binary value as values. In our specific 
            situation, the leaf names are usually UniProt IDs (i.e. A0A4P7DN87) 
            and the binary value either 1, 0 or -1 (nothing).

        Returns
        -------
        content : str
            Binary annotation file content.
        '''

        content = 'DATASET_BINARY\n'
        content += f'SEPARATOR {self.separator}\n'
        content += f'DATASET_LABEL{self.sep}binary\n'
        content += f'COLOR{self.sep}{self.color}\n'
        content += f'FIELD_LABELS{self.sep}{self.sep.join(self.field_labels)}\n'
        content += f'FIELD_COLORS{self.sep}{self.sep.join(self.field_colors)}\n'
        content += f'FIELD_SHAPES{self.sep}{self.sep.join(map(str, self.field_shapes))}\n'
        content += f'SHOW_INTERNAL{self.sep}{self.show_internal}\n'
        content += f'MARGIN{self.sep}{self.margin}\n'
        content += f'HEIGHT_FACTOR{self.sep}{self.height_factor}\n'
        content += f'SYMBOL_SPACING{self.sep}{self.symbol_spacing}\n'
        content += 'DATA\n'
        for branch in data:
            binary = self.sep.join(map(str, data[branch]))
            content += f'{branch}{self.sep}{binary}\n'
        
        return content