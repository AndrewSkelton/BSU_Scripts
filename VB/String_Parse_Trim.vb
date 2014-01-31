Sub parse_chr_extract()
    Dim tmp_str_ As String
    'Dim output_str As String
    For Each c In Range("A2:A174").Cells
        tmp_str_ = c.Text
        Range("B" & c.Row).Value = Left(tmp_str_, 3)
        Range("C" & c.Row).Value = Right(Left(tmp_str_, 4), 1)
    Next
End Sub
