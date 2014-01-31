Sub parse_split()
    Dim tmp_str_ As String
    Dim str_split_() As String
    For Each c In Range("A2:A174").Cells
        tmp_str_ = c.Value
        str_split_ = Split(tmp_str_, " ")
        Range("B" & c.Row).Value = str_split_(1)
        Range("C" & c.Row).Value = str_split_(2)
    Next
End Sub
