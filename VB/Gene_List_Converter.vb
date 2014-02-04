Sub Gene_Lists2()
    Dim MyObj As Object, MySource As Object, file As Variant, workingDir As String
    Dim pointer_ As Integer
    pointer_ = 2
    workingDir = "C:\Users\Andrew\Documents\Work\Young\Williams\Diff_Exp_Genes\"
    file = Dir(workingDir)
    While (file <> "")
        If Right(file, 3) = "txt" Then
            Filename_ = Replace(file, ".txt", "")
            
            Workbooks.Add
            With ActiveSheet.QueryTables.Add(Connection:= _
                "TEXT;" & workingDir & file _
                , Destination:=Range("$A$1"))
                .Name = Filename_
                .FieldNames = True
                .RowNumbers = False
                .FillAdjacentFormulas = False
                .PreserveFormatting = True
                .RefreshOnFileOpen = False
                .RefreshStyle = xlInsertDeleteCells
                .SavePassword = False
                .SaveData = True
                .AdjustColumnWidth = True
                .RefreshPeriod = 0
                .TextFilePromptOnRefresh = False
                .TextFilePlatform = 850
                .TextFileStartRow = 1
                .TextFileParseType = xlDelimited
                .TextFileTextQualifier = xlTextQualifierDoubleQuote
                .TextFileConsecutiveDelimiter = False
                .TextFileTabDelimiter = True
                .TextFileSemicolonDelimiter = False
                .TextFileCommaDelimiter = False
                .TextFileSpaceDelimiter = False
                .TextFileColumnDataTypes = Array(1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
                .TextFileTrailingMinusNumbers = True
                .Refresh BackgroundQuery:=False
            End With
            
            
            Range("A1").Select
            Range(Selection, Selection.End(xlToRight)).Select
            Selection.Cut Destination:=Range("B1:K1")
            Columns("A:A").Select
            Selection.Delete Shift:=xlToLeft
            
            For Each QTable In ActiveSheet.QueryTables
                QTable.Delete
            Next
            
            Range("A1").Select
            Range(Selection, Selection.End(xlDown)).Select
            Range(Selection, Selection.End(xlToRight)).Select

            ActiveSheet.ListObjects.Add(xlSrcRange, Range("$A$1:$J$" & Application.CountA(Columns(1))), , xlYes).Name = _
                "Table1"
            Range("Table1[#All]").Select
            ActiveSheet.ListObjects("Table1").TableStyle = "TableStyleMedium16"
            
            Columns("A:O").EntireColumn.AutoFit
            
            'save and close
            'switch back to workbook
            ActiveWorkbook.SaveAs Filename:= _
            workingDir & Filename_ & ".xlsx", FileFormat:= _
                xlOpenXMLWorkbook, CreateBackup:=False
            ActiveWindow.Close
            
            file = Dir
        End If
    Wend
End Sub
