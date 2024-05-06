'Excel.inc - Include file for BIFF 2.1 specifications to write Excel files.
'
'Converted from VB source to PowerBasic, November 2001 & February 2002.
'Paul Squires (2001-2002) support@planetsquires.com
'
'

'Class file for writing Microsoft Excel BIFF 2.1 files.

'This class is intended for users who do not want to use the huge
'Jet or ADO providers if they only want to export their data to
'an Excel compatible file.

'Newer versions of Excel use the OLE Structure Storage methods
'which are quite complicated.

'Paul Squires, February 25, 2001
'support@planetsquires.com


'constants to hold cell alignment
    %xlsGeneralAlign = 0
    %xlsLeftAlign = 1
    %xlsCentreAlign = 2
    %xlsRightAlign = 3
    %xlsFillCell = 4
    %xlsLeftBorder = 8
    %xlsRightBorder = 16
    %xlsTopBorder = 32
    %xlsBottomBorder = 64
    %xlsShaded = 128

'constants to handle selecting the font for the cell
    'used by rgbAttr2
    'bits 0-5 handle the *picture* formatting, not bold/underline etc...
    'bits 6-7 handle the font number
    %xlsFont0 = 0
    %xlsFont1 = 64
    %xlsFont2 = 128
    %xlsFont3 = 192

    'used by rgbAttr1
    'bits 0-5 must be zero
    'bit 6 locked/unlocked
    'bit 7 hidden/not hidden
    %xlsCellNormal = 0
    %xlsCellLocked = 64
    %xlsCellHidden = 128


'set up variables to hold the spreadsheet's layout
    %xlsLeftMargin = 38
    %xlsRightMargin = 39
    %xlsTopMargin = 40
    %xlsBottomMargin = 41


   'add these enums together. For example: xlsBold + xlsUnderline
   %xlsNoFormat = 0
   %xlsBold = 1
   %xlsItalic = 2
   %xlsUnderline = 4
   %xlsStrikeout = 8


TYPE FONT_RECORD
   opcode AS INTEGER  '49
   length AS INTEGER  '5+len(fontname)
   FontHeight AS INTEGER

   'bit0 bold, bit1 italic, bit2 underline, bit3 strikeout, bit4-7 reserved
   FontAttributes1 AS BYTE
   FontAttributes2 AS BYTE  'reserved - always 0
   FontNameLength AS BYTE
END TYPE


TYPE PASSWORD_RECORD
   opcode AS INTEGER  '47
   length AS INTEGER  'len(password)
END TYPE


TYPE HEADER_FOOTER_RECORD
   opcode AS INTEGER  '20 Header, 21 Footer
   length AS INTEGER  '1+len(text)
   TextLength AS BYTE
END TYPE


TYPE PROTECT_SPREADSHEET_RECORD
   opcode AS INTEGER  '18
   length AS INTEGER  '2
   Protect AS INTEGER
END TYPE

TYPE FORMAT_COUNT_RECORD
   opcode AS INTEGER  '1f
   length AS INTEGER '2
   Count AS INTEGER
END TYPE

TYPE FORMAT_RECORD
   opcode AS INTEGER  '1e
   length AS INTEGER  '1+len(format)
   FormatLength AS BYTE 'len(format)
END TYPE '+ followed by the Format-Picture

TYPE COLWIDTH_RECORD
   opcode AS INTEGER  '36
   length AS INTEGER  '4
   col1 AS BYTE       'first column
   col2 AS BYTE       'last column
   ColumnWidth AS INTEGER   'at 1/256th of a character
END TYPE

'Beginning Of File record
TYPE BEG_FILE_RECORD
  opcode AS INTEGER
  length AS INTEGER
  version AS INTEGER
  ftype AS INTEGER
END TYPE

'End Of File record
TYPE END_FILE_RECORD
  opcode AS INTEGER
  length AS INTEGER
END TYPE

'true/false to print gridlines
TYPE PRINT_GRIDLINES_RECORD
  opcode AS INTEGER
  length AS INTEGER
  PrintFlag AS INTEGER
END TYPE

'Integer record
TYPE tInteger
  opcode AS INTEGER
  length AS INTEGER
  Row AS INTEGER     'unsigned integer
  col AS INTEGER

  'rgbAttr1 handles whether cell is hidden and/or locked
  rgbAttr1 AS BYTE

  'rgbAttr2 handles the Font# and Formatting assigned to this cell
  rgbAttr2 AS BYTE

  'rgbAttr3 handles the Cell Alignment/borders/shading
  rgbAttr3 AS BYTE

  intValue AS INTEGER  'the actual integer value
END TYPE

'Number record
TYPE tNumber
  opcode AS INTEGER
  length AS INTEGER
  Row AS INTEGER
  col AS INTEGER
  rgbAttr1 AS BYTE
  rgbAttr2 AS BYTE
  rgbAttr3 AS BYTE
  NumberValue AS DOUBLE  '8 Bytes
END TYPE

'Label (Text) record
TYPE tText
  opcode AS INTEGER
  length AS INTEGER
  Row AS INTEGER
  col AS INTEGER
  rgbAttr1 AS BYTE
  rgbAttr2 AS BYTE
  rgbAttr3 AS BYTE
  TextLength AS BYTE
END TYPE

TYPE MARGIN_RECORD_LAYOUT
  opcode AS INTEGER
  length AS INTEGER
  MarginValue AS DOUBLE  '8 bytes
END TYPE

TYPE HPAGE_BREAK_RECORD
  opcode AS INTEGER
  length AS INTEGER
  NumPageBreaks AS INTEGER
END TYPE

TYPE DEF_ROWHEIGHT_RECORD
  opcode AS INTEGER
  length AS INTEGER
  RowHeight AS INTEGER
END TYPE

TYPE ROW_HEIGHT_RECORD
  opcode AS INTEGER  '08
  length AS INTEGER  'should always be 16 bytes
  RowNumber AS INTEGER
  FirstColumn AS INTEGER
  LastColumn AS INTEGER
  RowHeight AS INTEGER  'written to file as 1/20ths of a point
  internal AS INTEGER
  DefaultAttributes AS BYTE  'set to zero for no default attributes
  FileOffset AS INTEGER
  rgbAttr1 AS BYTE
  rgbAttr2 AS BYTE
  rgbAttr3 AS BYTE
END TYPE

'GLOBAL xlsFileNumber AS LONG
GLOBAL xlsBufferSize AS LONG  'if > 0 then buffer is active, also holds size of buffer.

'create an array that will hold the rows where a horizontal page
'break will be inserted just before.
GLOBAL xlsHorizPageBreakRows() AS LONG
GLOBAL xlsNumHorizPageBreaks AS LONG


DECLARE FUNCTION xlsCreateFile(mFileName$, xlsFileNumber&) AS LONG
DECLARE FUNCTION xlsCloseFile(xlsFileNumber&) AS LONG
DECLARE FUNCTION xlsInsertHorizPageBreak(BYVAL lrow AS LONG) AS LONG
DECLARE FUNCTION xlsWriteInteger(BYVAL value%, BYVAL lrow&, BYVAL lcol&, BYVAL CellFont&, BYVAL CellAlignment&, BYVAL HiddenLocked&, BYVAL CellFormat&, BYVAL xlsFileNumber&) AS LONG
DECLARE FUNCTION xlsWriteNumber(BYVAL value#, BYVAL lrow&, BYVAL lcol&, BYVAL CellFont&, BYVAL CellAlignment&, BYVAL HiddenLocked&, BYVAL CellFormat&, BYVAL xlsFileNumber&) AS LONG
DECLARE FUNCTION xlsWriteText(value$, BYVAL lrow&, BYVAL lcol&, BYVAL CellFont&, BYVAL CellAlignment&, BYVAL HiddenLocked&, BYVAL CellFormat&, BYVAL xlsFileNumber&) AS LONG
DECLARE FUNCTION xlsWriteDate(DateString$, BYVAL lrow&, BYVAL lcol&, BYVAL CellFont&, BYVAL CellAlignment&, BYVAL HiddenLocked&, BYVAL CellFormat&, BYVAL xlsFileNumber&) AS LONG
DECLARE FUNCTION xlsSetMargin(BYVAL Margin&, BYVAL MarginValue#, BYVAL xlsFileNumber&) AS LONG
DECLARE FUNCTION xlsSetColumnWidth(BYVAL FirstColumn&, BYVAL LastColumn&, BYVAL WidthValue&, BYVAL xlsFileNumber&) AS LONG
DECLARE FUNCTION xlsSetFont(FontName$, BYVAL FontHeight&, BYVAL FontFormat&, BYVAL xlsFileNumber&) AS LONG
DECLARE FUNCTION xlsSetHeader(HeaderText$, xlsFileNumber&) AS LONG
DECLARE FUNCTION xlsSetFooter(FooterText$, xlsFileNumber&) AS LONG
DECLARE FUNCTION xlsSetFilePassword(PasswordText$, xlsFileNumber&) AS LONG
DECLARE FUNCTION xlsPrintGridLines(BYVAL TrueFalse&, BYVAL xlsFileNumber&) AS LONG
DECLARE FUNCTION xlsProtectSpreadsheet(BYVAL TrueFalse&, BYVAL xlsFileNumber&) AS LONG
DECLARE FUNCTION xlsWriteDefaultFormats(xlsFileNumber&) AS LONG
DECLARE FUNCTION xlsSetDefaultRowHeight(BYVAL HeightValue&, BYVAL xlsFileNumber&) AS LONG
DECLARE FUNCTION xlsSetRowHeight(BYVAL lrow&, BYVAL HeightValue&, BYVAL xlsFileNumber&) AS LONG
DECLARE FUNCTION ConvertRow(BYVAL lrow AS LONG) AS INTEGER
DECLARE FUNCTION ConvertCol(BYVAL lcol AS LONG) AS INTEGER
DECLARE FUNCTION DateToJulian&(DateString$) AS LONG
DECLARE FUNCTION CTOD(PBDate AS STRING) AS STRING
DECLARE FUNCTION xlsBuffer(BYVAL TrueFalse&, BYVAL BufferSize&) AS LONG
DECLARE FUNCTION UpdateBuffer(BufferString$, xlsFileNumber&) AS LONG



FUNCTION xlsCreateFile(mFileName$, xlsFileNumber&) AS LONG

    IF DIR$(mFileName$) > "" THEN
       KILL mFileName$
       IF ERRCLEAR THEN
          FUNCTION = -1
          EXIT FUNCTION
       END IF
    END IF

    DIM BEG_FILE_MARKER AS BEG_FILE_RECORD
    'beginning of file
    BEG_FILE_MARKER.opcode = 9
    BEG_FILE_MARKER.length = 4
    BEG_FILE_MARKER.version = 2
    BEG_FILE_MARKER.ftype = 10

    xlsFileNumber& = FREEFILE
    OPEN mFileName$ FOR BINARY AS #xlsFileNumber&

    'if the buffer us active then save the data to the buffer
    'otherwise then simply write to the file.
    IF xlsBufferSize THEN
        stat& = UpdateBuffer((BEG_FILE_MARKER), xlsFileNumber&)
    ELSE
        PUT #xlsFileNumber&, , BEG_FILE_MARKER  'must always be written first
        IF ERRCLEAR THEN
           FUNCTION = -1
           EXIT FUNCTION
        END IF
    END IF

    'write the default formats to the file
    'and return if error occured.
    IF xlsWriteDefaultFormats(xlsFileNumber&) THEN EXIT FUNCTION

    'create the Horizontal Page Break array
    REDIM xlsHorizPageBreakRows(0)
    xlsNumHorizPageBreaks = 0

    FUNCTION = 0  'return with no error

END FUNCTION



FUNCTION xlsCloseFile(xlsFileNumber&) AS LONG

    IF xlsFileNumber& = 0 THEN
       FUNCTION = -1
       EXIT FUNCTION
    END IF

    'write the horizontal page breaks if necessary
    IF xlsNumHorizPageBreaks > 0 THEN
       'the Horizontal Page Break array must be in sorted order.
       'Use a simple Bubble sort because the size of this array would
       'be pretty small most of the time. A QuickSort would probably
       'be overkill.
         DIM lLoop1 AS LONG
         DIM lLoop2 AS LONG
         DIM lTemp AS LONG
         FOR lLoop1 = UBOUND(xlsHorizPageBreakRows) TO LBOUND(xlsHorizPageBreakRows) STEP -1
           FOR lLoop2 = LBOUND(xlsHorizPageBreakRows) + 1 TO lLoop1
             IF xlsHorizPageBreakRows(lLoop2 - 1) > xlsHorizPageBreakRows(lLoop2) THEN
               lTemp = xlsHorizPageBreakRows(lLoop2 - 1)
               xlsHorizPageBreakRows(lLoop2 - 1) = xlsHorizPageBreakRows(lLoop2)
               xlsHorizPageBreakRows(lLoop2) = lTemp
             END IF
           NEXT lLoop2
         NEXT lLoop1

       'write the Horizontal Page Break Record
        DIM HORIZ_PAGE_BREAK AS HPAGE_BREAK_RECORD
        HORIZ_PAGE_BREAK.opcode = 27
        HORIZ_PAGE_BREAK.length = 2 + (xlsNumHorizPageBreaks * 2)
        HORIZ_PAGE_BREAK.NumPageBreaks = xlsNumHorizPageBreaks

        IF xlsBufferSize THEN
            stat& = UpdateBuffer((HORIZ_PAGE_BREAK), BYVAL xlsFileNumber&)
        ELSE
            PUT #xlsFileNumber&, , HORIZ_PAGE_BREAK
            IF ERRCLEAR THEN
               FUNCTION = -1
               EXIT FUNCTION
            END IF
        END IF

        'now write the actual page break values
        FOR x& = 1 TO UBOUND(xlsHorizPageBreakRows)
           st$ = MKI$(xlsHorizPageBreakRows(x&))
           IF xlsBufferSize THEN
               stat& = UpdateBuffer(st$, xlsFileNumber&)
           ELSE
               PUT #xlsFileNumber&, , st$
               IF ERRCLEAR THEN
                  FUNCTION = -1
                  EXIT FUNCTION
               END IF
           END IF
        NEXT
    END IF

    DIM END_FILE_MARKER AS END_FILE_RECORD
    'end of file marker
    END_FILE_MARKER.opcode = 10

    IF xlsBufferSize THEN
        'set xlsBufferSize to -1 which will flag the UpdateBuffer routine
        'to flush the buffer.
        xlsBufferSize = -1
        stat& = UpdateBuffer("", xlsFileNumber&)
    END IF

    PUT #xlsFileNumber&, , END_FILE_MARKER

    CLOSE #xlsFileNumber&

    FUNCTION = 0  'return with no error code

END FUNCTION



FUNCTION xlsInsertHorizPageBreak(BYVAL lrow AS LONG) AS LONG

'the row and column values are written to the excel file as
'unsigned integers. Therefore, must convert the longs to integer.
    IF lrow > 32767 THEN
       Row% = CINT(lrow - 65536)
    ELSE
       Row% = CINT(lrow) - 1    'rows/cols in Excel binary file are zero based
    END IF

    xlsNumHorizPageBreaks = xlsNumHorizPageBreaks + 1
    REDIM PRESERVE xlsHorizPageBreakRows(xlsNumHorizPageBreaks)

    xlsHorizPageBreakRows(xlsNumHorizPageBreaks) = Row%

    FUNCTION = 0

END FUNCTION



FUNCTION ConvertRow(BYVAL lrow AS LONG) AS INTEGER
'the row and column values are written to the excel file as
'integers. Therefore, must convert the longs to integer.

    IF lrow > 32767 THEN
       FUNCTION = CINT(lrow - 65536)
    ELSE
       FUNCTION = CINT(lrow) - 1    'rows/cols in Excel binary file are zero based
    END IF

END FUNCTION


FUNCTION ConvertCol(BYVAL lcol AS LONG) AS INTEGER
'the row and column values are written to the excel file as
'integers. Therefore, must convert the longs to integer.

    IF lcol > 32767 THEN
       FUNCTION = CINT(lcol - 65536)
    ELSE
       FUNCTION = CINT(lcol) - 1    'rows/cols in Excel binary file are zero based
    END IF

END FUNCTION




FUNCTION xlsWriteInteger(BYVAL value%, BYVAL lrow&, BYVAL lcol&, BYVAL CellFont&, BYVAL CellAlignment&, BYVAL HiddenLocked&, BYVAL CellFormat&, BYVAL xlsFileNumber&) AS LONG

    'convert the row, col from LONG to INTEGER.
    Row% = ConvertRow(lrow&)
    Col% = ConvertCol(lcol&)

    DIM INTEGER_RECORD AS tInteger
    INTEGER_RECORD.opcode = 2
    INTEGER_RECORD.length = 9
    INTEGER_RECORD.Row = Row%
    INTEGER_RECORD.col = col%
    INTEGER_RECORD.rgbAttr1 = CBYT(HiddenLocked&)
    INTEGER_RECORD.rgbAttr2 = CBYT(CellFont& + CellFormat&)
    INTEGER_RECORD.rgbAttr3 = CBYT(CellAlignment&)
    INTEGER_RECORD.intValue = value%

    IF xlsBufferSize THEN
       stat& = UpdateBuffer((INTEGER_RECORD), xlsFileNumber&)
    ELSE
       PUT #xlsFileNumber&, , INTEGER_RECORD
       IF ERRCLEAR THEN
          FUNCTION = -1
          EXIT FUNCTION
       END IF
    END IF

    FUNCTION = 0   'return with no error

END FUNCTION



FUNCTION xlsWriteNumber(BYVAL value#, BYVAL lrow&, BYVAL lcol&, BYVAL CellFont&, BYVAL CellAlignment&, BYVAL HiddenLocked&, BYVAL CellFormat&, BYVAL xlsFileNumber&) AS LONG

    'convert the row, col from LONG to INTEGER.
    Row% = ConvertRow(lrow&)
    Col% = ConvertCol(lcol&)

    DIM NUMBER_RECORD AS tNumber
    NUMBER_RECORD.opcode = 3
    NUMBER_RECORD.length = 15
    NUMBER_RECORD.Row = Row%
    NUMBER_RECORD.col = col%
    NUMBER_RECORD.rgbAttr1 = CBYT(HiddenLocked&)
    NUMBER_RECORD.rgbAttr2 = CBYT(CellFont& + CellFormat&)
    NUMBER_RECORD.rgbAttr3 = CBYT(CellAlignment&)
    NUMBER_RECORD.NumberValue = value#

    IF xlsBufferSize THEN
       stat& = UpdateBuffer((NUMBER_RECORD), xlsFileNumber&)
    ELSE
       PUT #xlsFileNumber&, , NUMBER_RECORD
       IF ERRCLEAR THEN
          FUNCTION = -1
          EXIT FUNCTION
       END IF
    END IF

    FUNCTION = 0   'return with no error

END FUNCTION


FUNCTION xlsWriteText(value$, BYVAL lrow&, BYVAL lcol&, BYVAL CellFont&, BYVAL CellAlignment&, BYVAL HiddenLocked&, BYVAL CellFormat&, BYVAL xlsFileNumber&) AS LONG

    'convert the row, col from LONG to INTEGER.
    Row% = ConvertRow(lrow&)
    Col% = ConvertCol(lcol&)

    DIM b AS BYTE
    st$ = value$
    l& = LEN(st$)

    DIM TEXT_RECORD AS tText
    TEXT_RECORD.opcode = 4
    TEXT_RECORD.length = 10
    'Length of the text portion of the record
    TEXT_RECORD.TextLength = l&

    'Total length of the record
    TEXT_RECORD.length = 8 + l&

    TEXT_RECORD.Row = Row%
    TEXT_RECORD.col = col%

    TEXT_RECORD.rgbAttr1 = CBYT(HiddenLocked&)
    TEXT_RECORD.rgbAttr2 = CBYT(CellFont& + CellFormat&)
    TEXT_RECORD.rgbAttr3 = CBYT(CellAlignment&)

    'Put record header
    IF xlsBufferSize THEN
       stat& = UpdateBuffer((TEXT_RECORD), xlsFileNumber&)
    ELSE
       PUT #xlsFileNumber&, , TEXT_RECORD
    END IF

    'Then the actual string data
    FOR a& = 1 TO l&
        b = ASC(MID$(st$, a&, 1))
        IF xlsBufferSize THEN
           stat& = UpdateBuffer(MKBYT$(b), xlsFileNumber&)
        ELSE
           PUT #xlsFileNumber&, , b
           IF ERRCLEAR THEN
              FUNCTION = -1
              EXIT FUNCTION
           END IF
        END IF
    NEXT


    FUNCTION = 0   'return with no error

END FUNCTION


FUNCTION xlsWriteDate(DateString$, BYVAL lrow&, BYVAL lcol&, BYVAL CellFont&, BYVAL CellAlignment&, BYVAL HiddenLocked&, BYVAL CellFormat&, BYVAL xlsFileNumber&) AS LONG

    'convert the row, col from LONG to INTEGER.
    Row% = ConvertRow(lrow&)
    Col% = ConvertCol(lcol&)

    'convert the DateString$ from YYYYMMDD to a Julian date number
    value&= (DateToJulian&(DateString$) - DateToJulian&("19000100")) + 1

    DIM NUMBER_RECORD AS tNumber
    NUMBER_RECORD.opcode = 3
    NUMBER_RECORD.length = 15
    NUMBER_RECORD.Row = Row%
    NUMBER_RECORD.col = col%
    NUMBER_RECORD.rgbAttr1 = CBYT(HiddenLocked&)
    NUMBER_RECORD.rgbAttr2 = CBYT(CellFont& + CellFormat&)
    NUMBER_RECORD.rgbAttr3 = CBYT(CellAlignment&)
    NUMBER_RECORD.NumberValue = CDBL(value&)

    IF xlsBufferSize THEN
        stat& = UpdateBuffer((NUMBER_RECORD), xlsFileNumber&)
    ELSE
        PUT #xlsFileNumber&, , NUMBER_RECORD
        IF ERRCLEAR THEN
           FUNCTION = -1
           EXIT FUNCTION
        END IF
    END IF

    FUNCTION = 0   'return with no error

END FUNCTION


FUNCTION xlsSetMargin(BYVAL Margin&, BYVAL MarginValue#, BYVAL xlsFileNumber&) AS LONG

    'write the spreadsheet's layout information (in inches)
    DIM MARGINRECORD AS MARGIN_RECORD_LAYOUT

    'Margin& should be one of the following....
    '%xlsLeftMargin = 38
    '%xlsRightMargin = 39
    '%xlsTopMargin = 40
    '%xlsBottomMargin = 41

    MARGINRECORD.opcode = Margin&
    MARGINRECORD.length = 8
    MARGINRECORD.MarginValue = MarginValue# 'in inches

    IF xlsBufferSize THEN
        stat& = UpdateBuffer((MARGINRECORD), xlsFileNumber&)
    ELSE
        PUT #xlsFileNumber&, , MARGINRECORD
        IF ERRCLEAR THEN
           FUNCTION = -1
           EXIT FUNCTION
        END IF
    END IF

    FUNCTION = 0

END FUNCTION


FUNCTION xlsSetColumnWidth(BYVAL FirstColumn&, BYVAL LastColumn&, BYVAL WidthValue&, BYVAL xlsFileNumber&) AS LONG

    DIM COLWIDTH AS COLWIDTH_RECORD
    COLWIDTH.opcode = 36
    COLWIDTH.length = 4
    COLWIDTH.col1 = CBYT(FirstColumn&) - 1
    COLWIDTH.col2 = CBYT(LastColumn&) - 1
    COLWIDTH.ColumnWidth = WidthValue& * 256  'values are specified as 1/256 of a character

    IF xlsBufferSize THEN
        stat& = UpdateBuffer((COLWIDTH), xlsFileNumber&)
    ELSE
        PUT #xlsFileNumber&, , COLWIDTH
        IF ERRCLEAR THEN
           FUNCTION = -1
           EXIT FUNCTION
        END IF
    END IF

    FUNCTION = 0

END FUNCTION



FUNCTION xlsSetFont(FontName$, BYVAL FontHeight&, BYVAL FontFormat&, BYVAL xlsFileNumber&) AS LONG

    'you can set up to 4 fonts in the spreadsheet file. When writing a value such
    'as a Text or Number you can specify one of the 4 fonts (numbered 0 to 3)

    DIM FONTNAME_RECORD AS FONT_RECORD

    l& = LEN(FontName$)

    FONTNAME_RECORD.opcode = 49
    FONTNAME_RECORD.length = 5 + l&
    FONTNAME_RECORD.FontHeight = FontHeight& * 20
    FONTNAME_RECORD.FontAttributes1 = CBYT(FontFormat&)  'bold/underline etc...
    FONTNAME_RECORD.FontAttributes2 = CBYT(0)            'reserved-always zero!!
    FONTNAME_RECORD.FontNameLength = CBYT(l&)

    IF xlsBufferSize THEN
        stat& = UpdateBuffer((FONTNAME_RECORD), xlsFileNumber&)
    ELSE
        PUT #xlsFileNumber&, , FONTNAME_RECORD
    END IF

    'Then the actual font name data
    DIM b AS BYTE
    FOR a& = 1 TO l&
       b = ASC(MID$(FontName$, a&, 1))
       IF xlsBufferSize THEN
           stat& = UpdateBuffer(MKBYT$(b), xlsFileNumber&)
       ELSE
           PUT #xlsFileNumber&, , b
           IF ERRCLEAR THEN
              FUNCTION = -1
              EXIT FUNCTION
           END IF
       END IF
    NEXT

    FUNCTION = 0

END FUNCTION



FUNCTION xlsSetHeader(HeaderText$, xlsFileNumber&) AS LONG

    DIM HEADER_RECORD AS HEADER_FOOTER_RECORD

    l& = LEN(HeaderText$)

    HEADER_RECORD.opcode = 20
    HEADER_RECORD.length = 1 + l&
    HEADER_RECORD.TextLength = CBYT(l&)

    IF xlsBufferSize THEN
       stat& = UpdateBuffer((HEADER_RECORD), xlsFileNumber&)
    ELSE
       PUT #xlsFileNumber&, , HEADER_RECORD
    END IF

    'Then the actual Header text
    DIM b AS BYTE
    FOR a& = 1 TO l&
       b = ASC(MID$(HeaderText$, a&, 1))
       IF xlsBufferSize THEN
          stat& = UpdateBuffer(MKBYT$(b), xlsFileNumber&)
       ELSE
          PUT #xlsFileNumber&, , b
          IF ERRCLEAR THEN
             FUNCTION = -1
             EXIT FUNCTION
          END IF
       END IF
    NEXT


    FUNCTION = 0

END FUNCTION



FUNCTION xlsSetFooter(FooterText$, xlsFileNumber&) AS LONG

    DIM FOOTER_RECORD AS HEADER_FOOTER_RECORD

    l& = LEN(FooterText$)

    FOOTER_RECORD.opcode = 21
    FOOTER_RECORD.length = 1 + l&
    FOOTER_RECORD.TextLength = CBYT(l&)

    IF xlsBufferSize THEN
       stat& = UpdateBuffer((FOOTER_RECORD), xlsFileNumber&)
    ELSE
       PUT #xlsFileNumber&, , FOOTER_RECORD
    END IF

    'Then the actual Header text
    DIM b AS BYTE
    FOR a& = 1 TO l&
       b = ASC(MID$(FooterText$, a&, 1))
       IF xlsBufferSize THEN
          stat& = UpdateBuffer(MKBYT$(b), xlsFileNumber&)
       ELSE
          PUT #xlsFileNumber&, , b
          IF ERRCLEAR THEN
             FUNCTION = -1
             EXIT FUNCTION
          END IF
       END IF
    NEXT

    FUNCTION = 0

END FUNCTION



FUNCTION xlsSetFilePassword(PasswordText$, xlsFileNumber&) AS LONG

    DIM FILE_PASSWORD_RECORD AS PASSWORD_RECORD

    l& = LEN(PasswordText$)

    FILE_PASSWORD_RECORD.opcode = 47
    FILE_PASSWORD_RECORD.length = l&

    IF xlsBufferSize THEN
       stat& = UpdateBuffer((FILE_PASSWORD_RECORD), xlsFileNumber&)
    ELSE
       PUT #xlsFileNumber&, , FILE_PASSWORD_RECORD
    END IF

    'Then the actual Password text
    DIM b AS BYTE
    FOR a& = 1 TO l&
       b = ASC(MID$(PasswordText$, a&, 1))
       IF xlsBufferSize THEN
          stat& = UpdateBuffer(MKBYT$(b), xlsFileNumber&)
       ELSE
          PUT #xlsFileNumber&, , b
          IF ERRCLEAR THEN
             FUNCTION = -1
             EXIT FUNCTION
          END IF
       END IF
    NEXT

    FUNCTION = 0

END FUNCTION




FUNCTION xlsPrintGridLines(BYVAL TrueFalse&, BYVAL xlsFileNumber&) AS LONG

    DIM GRIDLINES_RECORD AS PRINT_GRIDLINES_RECORD

    GRIDLINES_RECORD.opcode = 43
    GRIDLINES_RECORD.length = 2

    IF TrueFalse& = 0 THEN
       GRIDLINES_RECORD.PrintFlag = 0
    ELSE
       GRIDLINES_RECORD.PrintFlag = 1
    END IF

    IF xlsBufferSize THEN
       stat& = UpdateBuffer((GRIDLINES_RECORD), xlsFileNumber&)
    ELSE
       PUT #xlsFileNumber&, , GRIDLINES_RECORD

       IF ERRCLEAR THEN
          FUNCTION = -1
          EXIT FUNCTION
       END IF
    END IF

    FUNCTION = 0

END FUNCTION




FUNCTION xlsProtectSpreadsheet(BYVAL TrueFalse&, BYVAL xlsFileNumber&) AS LONG

    DIM PROTECT_RECORD AS PROTECT_SPREADSHEET_RECORD

    PROTECT_RECORD.opcode = 18
    PROTECT_RECORD.length = 2

    IF TrueFalse& = 0 THEN
       PROTECT_RECORD.Protect = 0
    ELSE
       PROTECT_RECORD.Protect = 1
    END IF

    IF xlsBufferSize THEN
       stat& = UpdateBuffer((PROTECT_RECORD), xlsFileNumber&)
    ELSE
       PUT #xlsFileNumber&, , PROTECT_RECORD

       IF ERRCLEAR THEN
          FUNCTION = -1
          EXIT FUNCTION
       END IF
    END IF

    FUNCTION = 0

END FUNCTION


FUNCTION xlsWriteDefaultFormats(xlsFileNumber&) AS LONG

    DIM cFORMAT_COUNT_RECORD AS FORMAT_COUNT_RECORD
    DIM cFORMAT_RECORD AS FORMAT_RECORD
    DIM lIndex AS LONG
    DIM aFormat(0 TO 23) AS STRING
    DIM l AS LONG
    DIM q AS STRING
    q = CHR$(34)

    aFormat(0) = "General"
    aFormat(1) = "0"
    aFormat(2) = "0.000"
    aFormat(3) = "#,##0"
    aFormat(4) = "#,##0.00"
    aFormat(5) = "#,##0\ " & q & "$" & q & ";\-#,##0\ " & q & "$" & q
    aFormat(6) = "#,##0\ " & q & "$" & q & ";[Red]\-#,##0\ " & q & "$" & q
    aFormat(7) = "#,##0.00\ " & q & "$" & q & ";\-#,##0.00\ " & q & "$" & q
    aFormat(8) = "#,##0.00\ " & q & "$" & q & ";[Red]\-#,##0.00\ " & q & "$" & q
    aFormat(9) = "0%"
    aFormat(10) = "0.00%"
    aFormat(11) = "0.00E+00"
    aFormat(12) = "yyyy-mm-dd"
    aFormat(13) = "dd/\ mmm\ yy"
    aFormat(14) = "dd/\ mmm"
    aFormat(15) = "mmm\ yy"
    aFormat(16) = "h:mm\ AM/PM"
    aFormat(17) = "h:mm:ss\ AM/PM"
    aFormat(18) = "hh:mm"
    aFormat(19) = "hh:mm:ss"
    aFormat(20) = "dd/mm/yy\ hh:mm"
    aFormat(21) = "##0.0E+0"
    aFormat(22) = "mm:ss"
    aFormat(23) = "@"

    cFORMAT_COUNT_RECORD.opcode = &H1F
    cFORMAT_COUNT_RECORD.length = &H2
    cFORMAT_COUNT_RECORD.Count = CINT(UBOUND(aFormat))

    IF xlsBufferSize THEN
       stat& = UpdateBuffer((cFORMAT_COUNT_RECORD), xlsFileNumber&)
    ELSE
       PUT #xlsFileNumber&, , cFORMAT_COUNT_RECORD
    END IF

    FOR lIndex = LBOUND(aFormat) TO UBOUND(aFormat)
        l = LEN(aFormat(lIndex))

        cFORMAT_RECORD.opcode = &H1E
        cFORMAT_RECORD.length = l + 1
        cFORMAT_RECORD.FormatLength = l

        IF xlsBufferSize THEN
           stat& = UpdateBuffer((cFORMAT_RECORD), xlsFileNumber&)
        ELSE
           PUT #xlsFileNumber&, , cFORMAT_RECORD
        END IF

        'Then the actual format
        DIM b AS BYTE, a AS LONG
        FOR a = 1 TO l
           b = ASC(MID$(aFormat(lIndex), a, 1))
           IF xlsBufferSize THEN
              stat& = UpdateBuffer(MKBYT$(b), xlsFileNumber&)
           ELSE
              PUT #xlsFileNumber&, , b
              IF ERRCLEAR THEN
                 FUNCTION = -1
                 EXIT FUNCTION
              END IF
           END IF
        NEXT
    NEXT

    FUNCTION = 0

END FUNCTION



FUNCTION xlsSetDefaultRowHeight(BYVAL HeightValue&, BYVAL xlsFileNumber&) AS LONG

'Height is defined in units of 1/20th of a point. Therefore, a 10-point font
'would be 200 (i.e. 200/20 = 10). This function takes a HeightValue such as
'14 point and converts it the correct size before writing it to the file.

    DIM DEFHEIGHT AS DEF_ROWHEIGHT_RECORD
    DEFHEIGHT.opcode = 37
    DEFHEIGHT.length = 2
    DEFHEIGHT.RowHeight = HeightValue& * 20  'convert points to 1/20ths of point

    IF xlsBufferSize THEN
       stat& = UpdateBuffer((DEFHEIGHT), xlsFileNumber&)
    ELSE
       PUT #xlsFileNumber&, , DEFHEIGHT

       IF ERRCLEAR THEN
          FUNCTION = -1
          EXIT FUNCTION
       END IF
    END IF

    FUNCTION = 0

END FUNCTION


FUNCTION xlsSetRowHeight(BYVAL lrow&, BYVAL HeightValue&, BYVAL xlsFileNumber&) AS LONG

    'convert the row, col from LONG to INTEGER.
    Row% = ConvertRow(lrow&)

    'Height is defined in units of 1/20th of a point. Therefore, a 10-point font
    'would be 200 (i.e. 200/20 = 10). This function takes a HeightValue such as
    '14 point and converts it the correct size before writing it to the file.

    DIM ROWHEIGHTREC AS ROW_HEIGHT_RECORD
    ROWHEIGHTREC.opcode = 8
    ROWHEIGHTREC.length = 16
    ROWHEIGHTREC.RowNumber = Row%
    ROWHEIGHTREC.FirstColumn = 0
    ROWHEIGHTREC.LastColumn = 256
    ROWHEIGHTREC.RowHeight = HeightValue& * 20 'convert points to 1/20ths of point
    ROWHEIGHTREC.internal = 0
    ROWHEIGHTREC.DefaultAttributes = 0
    ROWHEIGHTREC.FileOffset = 0
    ROWHEIGHTREC.rgbAttr1 = 0
    ROWHEIGHTREC.rgbAttr2 = 0
    ROWHEIGHTREC.rgbAttr3 = 0

    IF xlsBufferSize THEN
       stat& = UpdateBuffer((ROWHEIGHTREC), xlsFileNumber&)
    ELSE
       PUT #xlsFileNumber&, , ROWHEIGHTREC

       IF ERRCLEAR THEN
          FUNCTION = -1
          EXIT FUNCTION
       END IF
    END IF

    FUNCTION = 0

END FUNCTION


FUNCTION DateToJulian&(DateString$) AS LONG

  'DateString$ must be in YYYYMMDD

  LOCAL Elapsed AS LONG

  IF LEN(DateString$) <> 8 THEN
     FUNCTION = 0
     EXIT FUNCTION
  END IF

  Year& = VAL(LEFT$(DateString$, 4))
  month& = VAL(MID$(DateString$, 5, 2))
  day& = VAL(RIGHT$(DateString$, 2))

  IF month& < 3 THEN                      ' January or February?
     month& = month& + 12                   ' 13th or 14th month ....
     DECR year&                            ' .... of prev. year
  END IF

  Elapsed = INT((year& + 4712) * 365.25)  ' years elapsed
  Elapsed = Elapsed - (year& \ 100)       ' substract century leapdays
  Elapsed = Elapsed + (year& \ 400)       ' re-add valid ones
  Elapsed = Elapsed + _
            INT(30.6 * (month& - 1) + .2) ' months elapsed + adjustm.
  FUNCTION = Elapsed + day&               ' days of final month

END FUNCTION



FUNCTION CTOD$(PBDate AS STRING)
'¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
'Converts a PowerBASIC formated date string into an xBase date
' Usage:   ans$ = CTOD("12-13-2001")
'          ans$ = "20011213"
'¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
'First parse PBDate
 DIM month$, day$, year$
 DIM ChekYear&

  month$ = PARSE$(PBDate, ANY "-", 1)     'month
  day$ = PARSE$(PBDate, ANY "-", 2)       'day
  year$ = PARSE$(PBDate, ANY "-", 3)      'year

  ChekYear& = VAL(year$)

' Adjust if user enters only 2 digit year
   IF ChekYear& < 80   THEN     'assume we've rolled into the next century
      year$ = "20" + year$
   ELSEIF ChekYear& < 100 THEN   'between 1980 and 1999
      year$ = "19" + year$
   END IF

' Force the reformat of strings into proper xBase length.
  month$ = FORMAT$(VAL(month$), "00")
  day$ = FORMAT$(VAL(day$), "00")
  year$ = FORMAT$(VAL(year$), "0000")


  FUNCTION = (year$ & month$ & day$)

END FUNCTION


FUNCTION xlsBuffer(BYVAL TrueFalse&, BYVAL BufferSize&) AS LONG

    IF BufferSize& <= 0 THEN EXIT FUNCTION

    IF TrueFalse& <> 0 THEN
       xlsBufferSize = BufferSize&  'flag that buffer is active & also contains the buffer size
    END IF

END FUNCTION


FUNCTION UpdateBuffer(BufferData$, xlsFileNumber&) AS LONG

    STATIC FirstTime&
    STATIC Buffer$
    STATIC CurBufferPos&

    IF FirstTime& = 0 THEN
       'initilize the buffer
       Buffer$ = SPACE$(xlsBufferSize)
       CurBufferPos& = 0
       FirstTime& = -1
    END IF

    'if xlsBufferSize equals -1 then the flag to flush the buffer has
    'been set and the file is about to be closed.
    IF xlsBufferSize = -1 THEN
       PUT$ xlsFileNumber&, LEFT$(Buffer$, CurBufferPos&)
       FirstTime& = 0
       EXIT FUNCTION
    END IF

    'if the current data would be too big for the buffer then flush
    'the buffer and reset its size.
    IF LEN(BufferData$) + CurBufferPos& > xlsBufferSize& THEN
       PUT$ xlsFileNumber&, LEFT$(Buffer$, CurBufferPos&)
       CurBufferPos& = 0
    END IF

    'write the data to the buffer and update the buffer size indicator
    MID$(Buffer$, CurBufferPos& + 1) = BufferData$
    CurBufferPos& = CurBufferPos& + LEN(BufferData$)


END FUNCTION
