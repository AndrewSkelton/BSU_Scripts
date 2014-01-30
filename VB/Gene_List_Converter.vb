Sub Gene_Lists()
    
	Dim MyObj As Object, MySource As Object, file As Variant, workingDir As String
	workingDir = "C:\Users\Andrew\Documents\Work\Young\Shereen\Gene_Lists\txt\"
	file = Dir(workingDir)
	While (file <> "")
		If InStr(file, ".txt") > 0 Then
			file = Dir
			Filename_ = Replace(file, ".txt", "")
			MsgBox Filename_
			
			'Open new file and switch to it
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
				.TextFileColumnDataTypes = Array(1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
				.TextFileTrailingMinusNumbers = True
				.Refresh BackgroundQuery:=False
			End With
			
			Range("A1:O1").Select
			Selection.Cut Destination:=Range("B1:P1")
			
			Columns("A:A").Select
			Selection.Delete Shift:=xlToLeft
			Range("A1").Select
			Range(Selection, Selection.End(xlDown)).Select
			Range(Selection, Selection.End(xlToRight)).Select
			
			ActiveSheet.QueryTables(Filename_).Delete
			ActiveSheet.ListObjects.Add(xlSrcRange, Range("$A$1:$O$2"), , xlYes).Name = _
			"Table1"
			Range("Table1[#All]").Select
			ActiveSheet.ListObjects("Table1").TableStyle = "TableStyleMedium16"
			Columns("A:O").Select
			Range("Table1[[#Headers],[B]]").Activate
			Columns("A:O").EntireColumn.AutoFit
			
			'save and close
			'switch back to workbook
			 ActiveWorkbook.SaveAs Filename:= _
			workingDir & Filename_ & ".xlsx", FileFormat:= _
				xlOpenXMLWorkbook, CreateBackup:=False
			ActiveWindow.Close
		End If
	Wend
End Sub
