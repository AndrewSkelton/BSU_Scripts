sub stringBuilder()
    Dim tmp_str_ As String
    For Each c In Range("H2:H174").Cells
        tmp_str_ = tmp_str_ & """" & c.Value & "_" & Range("I" & c.Row).Value & _
                    "_" & Range("J" & c.Row).Value & "_" & Range("K" & c.Row).Value & _
                    """, "
    Next
    Range("M1").Value = tmp_str_
End Sub
