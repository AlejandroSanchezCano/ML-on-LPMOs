class Annotation():
    def __init__(
            self,
            separator = 'SPACE',
            name = 'color_strip1',
            color = '#000',
            margin = 0,
            show_internal = 0,

            color_branches = 0,
            strip_width = 70,
            border_width = 10,
            border_color = '#000',

            field_labels = ('C1', 'C4', 'substrate'),
            field_colors = ('#000', '#000', '#66f'),
            field_shapes = (5, 4, 3),
            height_factor = 85,
            symbol_spacing = 10

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

    def colorstrip(self, data):
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

    def binary(self, data):
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