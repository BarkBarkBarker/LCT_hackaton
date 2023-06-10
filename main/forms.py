from django import forms
from decimal import Decimal


class Form1(forms.Form):
    D1 = forms.DecimalField(label="", required=False,
                    widget=forms.TextInput(attrs={'placeholder': 'd'}))
    X1 = forms.DecimalField(label="", required=False,
                    widget=forms.TextInput(attrs={'placeholder': 'x'}))
    Y1 = forms.DecimalField(label="", required=False,
                    widget=forms.TextInput(attrs={'placeholder': 'y'}))
    Z1 = forms.DecimalField(label="", required=False,
                    widget=forms.TextInput(attrs={'placeholder': 'z'}))
    def __init__(self, *args, **kwargs):

        super(Form1, self).__init__(*args, **kwargs)
        for visible in self.visible_fields():
            visible.field.widget.attrs['class'] = 'form'


class Form2(forms.Form):
    D2 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'd'}))
    X2 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'x'}))
    Y2 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'y'}))
    Z2 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'z'}))

    def __init__(self, *args, **kwargs):
        super(Form2, self).__init__(*args, **kwargs)
        for visible in self.visible_fields():
            visible.field.widget.attrs['class'] = 'form'
            self.field_order = [
                'X2',
                'Y2',
                'Z2',
                'D2']
class Form3(forms.Form):
    D3 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'd'}))
    X3 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'x'}))
    Y3 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'y'}))
    Z3 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'z'}))

    def __init__(self, *args, **kwargs):
        super(Form3, self).__init__(*args, **kwargs)
        for visible in self.visible_fields():
            visible.field.widget.attrs['class'] = 'form'
            self.field_order = [
                'X3',
                'Y3',
                'Z3',
                'D3']
class Form4(forms.Form):
    Y21 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'y2'}))
    X21 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'x2'}))
    Y11 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'y1'}))
    X11 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'x1'}))
    def __init__(self, *args, **kwargs):
        super(Form4, self).__init__(*args, **kwargs)
        for visible in self.visible_fields():
            visible.field.widget.attrs['class'] = 'form'
        self.field_order = [
            'X11',
            'Y11',
            'X21',
            'Y21']

class Form5(forms.Form):
    Y22 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'y4'}))
    X22 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'x4'}))
    Y12 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'y3'}))
    X12 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'x3'}))
    def __init__(self, *args, **kwargs):
        super(Form5, self).__init__(*args, **kwargs)
        for visible in self.visible_fields():
            visible.field.widget.attrs['class'] = 'form'
            self.field_order = [
                'X12',
                'Y12',
                'X22',
                'Y22']

class Form6(forms.Form):
    Y23 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'y6'}))
    X23 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'x6'}))
    Y13 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'y5'}))
    X13 = forms.DecimalField(label="", required=False,
                            widget=forms.TextInput(attrs={'placeholder': 'x5'}))

    def __init__(self, *args, **kwargs):
        super(Form6, self).__init__(*args, **kwargs)
        for visible in self.visible_fields():
            visible.field.widget.attrs['class'] = 'form'
            self.field_order = [
                'X13',
                'Y13',
                'X23',
                'Y23']

class Form7(forms.Form):
    DW2 = forms.DecimalField(label="", required=False,
                             widget=forms.TextInput(attrs={'placeholder': 'D23'}))
    DW3 = forms.DecimalField(label="", required=False,
                             widget=forms.TextInput(attrs={'placeholder': 'D34'}))
    DW4 = forms.DecimalField(label="", required=False,
                             widget=forms.TextInput(attrs={'placeholder': 'D45'}))
    DW5 = forms.DecimalField(label="", required=False,
                             widget=forms.TextInput(attrs={'placeholder': 'D52'}))

    def __init__(self, *args, **kwargs):
        super(Form7, self).__init__(*args, **kwargs)
        for visible in self.visible_fields():
            visible.field.widget.attrs['class'] = 'form'
            self.field_order = [
                'DW2',
                'DW3',
                'DW4',
                'DW5']


class Form8(forms.Form):
    D4 = forms.DecimalField(label="", required=False)
    X4 = forms.DecimalField(label="", required=False)
    Y4 = forms.DecimalField(label="", required=False)
    Z4 = forms.DecimalField(label="", required=False)

    D5 = forms.DecimalField(label="", required=False)
    X5 = forms.DecimalField(label="", required=False)
    Y5 = forms.DecimalField(label="", required=False)
    Z5 = forms.DecimalField(label="", required=False)

    D6 = forms.DecimalField(label="", required=False)
    X6 = forms.DecimalField(label="", required=False)
    Y6 = forms.DecimalField(label="", required=False)
    Z6 = forms.DecimalField(label="", required=False)

    D7 = forms.DecimalField(label="", required=False)
    X7 = forms.DecimalField(label="", required=False)
    Y7 = forms.DecimalField(label="", required=False)
    Z7 = forms.DecimalField(label="", required=False)

    D8 = forms.DecimalField(label="", required=False)
    X8 = forms.DecimalField(label="", required=False)
    Y8 = forms.DecimalField(label="", required=False)
    Z8 = forms.DecimalField(label="", required=False)

    def __init__(self, *args, **kwargs):
        super(Form8, self).__init__(*args, **kwargs)
        for visible in self.visible_fields():
            visible.field.widget.attrs['class'] = 'form'


class Form9(forms.Form):

    D9 = forms.DecimalField(label="", required=False)
    X9 = forms.DecimalField(label="", required=False)
    Y9 = forms.DecimalField(label="", required=False)
    Z9 = forms.DecimalField(label="", required=False)

    D10 = forms.DecimalField(label="", required=False)
    X10 = forms.DecimalField(label="", required=False)
    Y10 = forms.DecimalField(label="", required=False)
    Z10 = forms.DecimalField(label="", required=False)

    D11 = forms.DecimalField(label="", required=False)
    X11 = forms.DecimalField(label="", required=False)
    Y11 = forms.DecimalField(label="", required=False)
    Z11 = forms.DecimalField(label="", required=False)

    D12 = forms.DecimalField(label="", required=False)
    X12 = forms.DecimalField(label="", required=False)
    Y12 = forms.DecimalField(label="", required=False)
    Z12 = forms.DecimalField(label="", required=False)

    def __init__(self, *args, **kwargs):
        super(Form9, self).__init__(*args, **kwargs)
        for visible in self.visible_fields():
            visible.field.widget.attrs['class'] = 'form'


class HelpfulForm(forms.Form):
    form_template_name = "/templates/main/form_snippet.html"
    check_enable = forms.BooleanField(required=False)
    num = forms.DecimalField(required=False)

    def __init__(self, *args, **kwargs):
        super(HelpfulForm, self).__init__(*args, **kwargs)
        for visible in self.visible_fields():
            visible.field.widget = visible.field.hidden_widget()
