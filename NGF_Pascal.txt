1. Initialization
procedure TFrmMain.Btn_InitClick(Sender: TObject);
var
  i, j, k, Begin_pos, End_Pos, IPos, index, totalnum, Samples, L, m, a, b, a_b: integer;
  neighbors, Domain_weighted, Diff1, Diff2, P1, P2, Complex_Weight, max_dom, min_dom: double;
  s, s1: string;
  LoadFile: TStringList;
  NodeInfo: array of NodeInformation;
  ChangeNode: NodeInformation;
  sub_score: array of double;
  Weighted: Boolean;
  Detect_Node: array of array of integer; //复合物
  Matrix_Complex: array of array of integer;
  Select_Proteins: array of integer;
  Neighbor_list: array of array of integer; //邻接节点集合
begin
  if FileExists(Edit2.Text) = False then
  begin
    showmessage('PPI network is not Found!');
    Exit;
  end;
  Btn_Init.Enabled := False;
  Btn_NGF.Enabled := False;
  //------------------------------------------------
  LoadFile := TStringList.Create;
  Nodes := TStringList.Create;
  Nodes.Clear;
  Nodes.Sorted := True;
  LoadFile.LoadFromFile(Edit2.Text);
  Gauge1.MaxValue := LoadFile.Count - 1;
  Gauge1.MinValue := 0;
  Gauge1.Progress := 0;
  for i := 0 to LoadFile.Count - 1 do
  begin
    Gauge1.Progress := Gauge1.Progress + 1;
    s := LoadFile.Strings[i];
    IPos := pos(#9, s);
    s1 := Copy(s, 1, IPos - 1);
    Begin_pos := Nodes.IndexOf(s1);
    if Begin_pos < 0 then
      Nodes.Add(s1);
    Delete(s, 1, IPos);
    IPos := pos(#9, s);
    if IPos = 0 then
    begin
      End_Pos := Nodes.IndexOf(s);
      s1 := s;
    end
    else
    begin
      s1 := Copy(s, 1, IPos - 1);
      End_Pos := Nodes.IndexOf(s1);
    end;
    if End_Pos < 0 then
      Nodes.Add(s1);
    Application.ProcessMessages;
  end;
  //-----------------------------------------------------------------
  setlength(Edge_Info, 0);
  setlength(Edge_Info, LoadFile.Count);
  setlength(MatrixEdge, 0);
  setlength(MatrixEdge, nodes.Count);
  setlength(Neighbor_list, nodes.Count);
  for i := 0 to nodes.Count - 1 do
  begin
    setlength(MatrixEdge[i], nodes.Count);
    setlength(Neighbor_list[i], 1);
    Neighbor_list[i][0] := i;
  end;
  setlength(NodeInfo, 0);
  setlength(NodeInfo, Nodes.Count);
  Weighted := False;
  for i := 0 to LoadFile.Count - 1 do
  begin
    s := LoadFile.Strings[i];
    IPos := pos(#9, s);
    s1 := Copy(s, 1, IPos - 1);
    Begin_pos := Nodes.IndexOf(s1);
    if NodeInfo[Begin_pos].Node = '' then
    begin
      NodeInfo[Begin_pos].Node := s1;
      NodeInfo[Begin_pos].degree := 0;
      NodeInfo[Begin_pos].GoNums := 0;
    end;
    Delete(s, 1, IPos);
    IPos := pos(#9, s);
    if IPos = 0 then
    begin
      End_Pos := Nodes.IndexOf(s);
      s1 := s;
    end
    else
    begin
      s1 := Copy(s, 1, IPos - 1);
      Delete(s, 1, IPos);
      End_Pos := Nodes.IndexOf(s1);
      Weighted := True;
    end;
    if NodeInfo[End_Pos].Node = '' then
    begin
      NodeInfo[End_Pos].Node := s1;
      NodeInfo[End_Pos].degree := 0;
      NodeInfo[End_Pos].GoNums := 0;
    end;
    Edge_Info[i].BeginNode := Begin_pos;
    Edge_Info[i].EndNode := End_Pos;
    if Weighted = True then
      Edge_Info[i].ECC := StrToFloat(s)
    else
      Edge_Info[i].ECC := 0;
    MatrixEdge[Begin_pos][End_Pos].Flag := 1;
    MatrixEdge[End_Pos][Begin_pos].Flag := 1;
    NodeInfo[End_Pos].degree := NodeInfo[End_Pos].degree + 1;
    NodeInfo[Begin_pos].degree := NodeInfo[Begin_pos].degree + 1;
    Application.ProcessMessages;
  end;
  Total_Proteins := length(NodeInfo);
    //-------------------------------------------------------------
  GO := TStringList.Create;
  GO.Clear;
  GO.Sorted := True;
  if Radiobutton1.Checked then
    LoadFile.LoadFromFile(ExtractFilePath(Application.ExeName) + 'GO_P.txt')
  else if Radiobutton2.Checked then
    LoadFile.LoadFromFile(ExtractFilePath(Application.ExeName) + 'GO_F.txt')
  else
    LoadFile.LoadFromFile(ExtractFilePath(Application.ExeName) + 'GO_C.txt');
  Gauge1.MaxValue := LoadFile.Count;
  Gauge1.MinValue := 0;
  Gauge1.Progress := 0;
  k := 0;
  for i := 0 to LoadFile.Count - 1 do
  begin
    Gauge1.Progress := Gauge1.Progress + 1;
    s := LoadFile.Strings[i];
    IPos := pos(#9, s);
    s1 := Copy(s, 1, IPos - 1);
    Begin_pos := Nodes.IndexOf(s1);
    if Begin_pos < 0 then
    begin
      continue;
    end;
    Delete(s, 1, IPos);
    End_Pos := GO.IndexOf(s);
    if End_Pos < 0 then
      GO.Add(s);
    Application.ProcessMessages;
  end;
  setlength(Matrix_GO, 0);
  setlength(Matrix_GO, Nodes.Count, GO.Count);
  Gauge1.MaxValue := LoadFile.Count;
  Gauge1.MinValue := 0;
  Gauge1.Progress := 0;
  Memo3.Lines.Clear;
  for i := 0 to LoadFile.Count - 1 do
  begin
    Gauge1.Progress := Gauge1.Progress + 1;
    s := LoadFile.Strings[i];
    IPos := pos(#9, s);
    s1 := Copy(s, 1, IPos - 1);
    Begin_pos := Nodes.IndexOf(s1);
    Delete(s, 1, IPos);
    End_Pos := GO.IndexOf(s);
    if (Begin_Pos >= 0) and (End_Pos >= 0) then
    begin
      Matrix_GO[Begin_pos][End_Pos] := 1;
      NodeInfo[Begin_pos].GONums := NodeInfo[Begin_pos].GONums + 1;
    end;
    Application.ProcessMessages;
  end;
  for i := 0 to nodes.Count - 1 do
    if nodeinfo[i].GoNums > 0 then
      Memo3.Lines.Add(inttostr(nodeinfo[i].GoNums));
  memo3.Lines.SaveToFile('GoNum.txt');
   //-------------------------------------------------------------
  setlength(Detect_Node, 0);
  LoadFile.LoadFromFile(ExtractFilePath(Application.ExeName) + 'Know_complexes.txt');
  for i := 0 to LoadFile.Count - 1 do
  begin
    s := LoadFile.Strings[i];
    if LeftStr(s, 7) = 'Complex' then
    begin
      s := ReverseString(s);
      j := 1;
      s1 := '';
      while true do
      begin
        if (s[j] = ' ') or (s[j] = #9) then
          break;
        s1 := s1 + s[j];
        Inc(j);
      end;
      s1 := ReverseString(s1);
      setlength(Detect_Node, length(Detect_Node) + 1);
      setlength(Detect_Node[high(Detect_Node)], StrToInt(s1));
      Index := 0;
    end
    else
    begin
      IPos := Nodes.IndexOf(s);
      if Nodes.IndexOf(s) >= 0 then
        Detect_Node[high(Detect_Node)][Index] := IPos;
      Inc(Index);
    end;
  end;

  setlength(Matrix_Complex, Nodes.Count, length(Detect_Node));
  for i := 0 to Nodes.Count - 1 do
  begin
    Index := 0;
    for j := 0 to high(Detect_Node) do
    begin
      for k := 0 to high(Detect_Node[j]) do
        if Detect_Node[j][k] = i then
          break;
      if k <= high(Detect_Node[j]) then
        Matrix_Complex[i][j] := 1
      else
        Matrix_Complex[i][j] := 0;

      Index := Index + Matrix_Complex[i][j];
    end;
    NodeInfo[i].ComplexNum := Index;
    Application.ProcessMessages;
  end;
   //--------------------------------------------------------------
  Domain := TStringList.Create;
  Domain.Clear;
  Domain.Sorted := True;
  LoadFile.LoadFromFile(ExtractFilePath(Application.ExeName) + 'Domain.txt');
  Gauge1.MaxValue := LoadFile.Count;
  Gauge1.MinValue := 0;
  Gauge1.Progress := 0;
  k := 0;
  for i := 0 to LoadFile.Count - 1 do
  begin
    Gauge1.Progress := Gauge1.Progress + 1;
    s := LoadFile.Strings[i];
    IPos := pos(#9, s);
    s1 := Copy(s, 1, IPos - 1);
    Begin_pos := Nodes.IndexOf(s1);
    if Begin_pos < 0 then
    begin
      continue;
    end;
    Delete(s, 1, IPos);
    End_Pos := Domain.IndexOf(s);
    if End_Pos < 0 then
      Domain.Add(s);
    Application.ProcessMessages;
  end;

  setlength(Domain_frequency, Domain.Count);
  for i := 0 to domain.Count - 1 do
    Domain_frequency[i] := 0;
  setlength(Matrix_Domain, 0);
  setlength(Matrix_Domain, Nodes.Count, Domain.Count);
  Gauge1.MaxValue := LoadFile.Count;
  Gauge1.MinValue := 0;
  Gauge1.Progress := 0;
  for i := 0 to LoadFile.Count - 1 do
  begin
    Gauge1.Progress := Gauge1.Progress + 1;
    s := LoadFile.Strings[i];
    IPos := pos(#9, s);
    s1 := Copy(s, 1, IPos - 1);
    Begin_pos := Nodes.IndexOf(s1);
    Delete(s, 1, IPos);
    End_Pos := Domain.IndexOf(s);
    if (Begin_Pos >= 0) and (End_Pos >= 0) then
    begin
      Matrix_Domain[Begin_pos][End_Pos] := 1;
      NodeInfo[Begin_pos].DomainNums := NodeInfo[Begin_pos].DomainNums + 1;
      inc(Domain_frequency[End_Pos]);
    end;
    Application.ProcessMessages;
  end;

  //--------------------------------------------
  Gauge1.MaxValue := high(Edge_Info);
  Gauge1.MinValue := 0;
  Gauge1.Progress := 0;

  for i := low(Edge_Info) to high(Edge_Info) do
  begin
    Gauge1.Progress := Gauge1.Progress + 1;
    //------------------------------ECC-----------------------------------------
    neighbors := 0; Diff1 := 0; Diff2 := 0;
    for k := low(NodeInfo) to high(NodeInfo) do
    begin
      if (k <> Edge_Info[i].BeginNode) and (k <> Edge_Info[i].EndNode) then
      begin
        if (MatrixEdge[Edge_Info[i].EndNode][k].Flag = 1) and (MatrixEdge[Edge_Info[i].BeginNode][k].Flag = 1) then
          neighbors := neighbors + 1
        else if (MatrixEdge[Edge_Info[i].EndNode][k].Flag = 1) and (MatrixEdge[Edge_Info[i].BeginNode][k].Flag = 0) then
          Diff1 := Diff1 + 1
        else if (MatrixEdge[Edge_Info[i].EndNode][k].Flag = 0) and (MatrixEdge[Edge_Info[i].BeginNode][k].Flag = 1) then
          Diff2 := Diff2 + 1;
      end;
      Application.ProcessMessages;
    end;

    if (Min(NodeInfo[Edge_Info[i].BeginNode].degree, NodeInfo[Edge_Info[i].EndNode].degree) > 1) and (Weighted = False) then
    begin
      Edge_Info[i].ECC := 4 * power(neighbors, 2) / ((NodeInfo[Edge_Info[i].BeginNode].degree + neighbors) * (NodeInfo[Edge_Info[i].EndNode].degree + neighbors));
    end;
    //-------------------------------------------------------------------------------
    Application.ProcessMessages;
  end;



  //----------------------------------------------------------------------------
  Memo2.Lines.Clear;
  setlength(Node_Info, Nodes.Count);
  //----------------------------------------------------------------------------

  setlength(Matrix_Edge, Nodes.Count);
  IPos := 0;
  Gauge1.MaxValue := high(NodeInfo);
  Gauge1.MinValue := 0;
  Gauge1.Progress := 0;
  for i := 0 to high(NodeInfo) do
  begin
    Gauge1.Progress := Gauge1.Progress + 1;
    begin
      Node_Info[IPos] := NodeInfo[i];
      setlength(Matrix_Edge[IPos], Nodes.Count);
      for j := 0 to Nodes.Count - 1 do
      begin
        Matrix_Edge[IPos][j][0] := 0;
        Matrix_Edge[IPos][j][1] := 0;
        Matrix_Edge[IPos][j][2] := 0;
        Application.ProcessMessages;
      end;
      inc(IPos);
    end;
    Application.ProcessMessages;
  end;


  Gauge1.MaxValue := high(Edge_Info);
  Gauge1.MinValue := 0;
  Gauge1.Progress := 0;
  Memo1.Lines.Clear;
  for i := low(Edge_Info) to high(Edge_Info) do
  begin
    Gauge1.Progress := Gauge1.Progress + 1;
    Begin_pos := Nodes.IndexOf(NodeInfo[Edge_Info[i].BeginNode].Node);
    End_pos := Nodes.IndexOf(NodeInfo[Edge_Info[i].EndNode].Node);

    Matrix_Edge[Begin_pos][End_Pos][0] := Edge_Info[i].ECC;
    Matrix_Edge[End_Pos][Begin_pos][0] := Edge_Info[i].ECC;

    Application.ProcessMessages;
  end;

  Gauge1.MaxValue := nodes.Count;
  Gauge1.MinValue := 0;
  Gauge1.Progress := 0;
  memo2.Lines.Clear;
  max_dom := 0;
  min_dom := 10000;
  for i := 0 to nodes.Count - 1 do
  begin
    Gauge1.Progress := Gauge1.Progress + 1;

    for j := i + 1 to nodes.Count - 1 do
    begin

      a := NodeInfo[i].DomainNums;
      b := NodeInfo[j].DomainNums; ;
      a_b := 0;
      for k := 0 to domain.Count - 1 do
        a_b := a_b + Matrix_Domain[i][k] * Matrix_Domain[j][k]; ;
      if a_b = 0 then
        Domain_weighted := 0
      else
      begin
        P1 := 1.0 * power(domain.Count, a) / (power(domain.Count, a_b) * power(domain.Count - a_b, a - a_b));
        p1 := log10(p1);
        p2 := 1.0 * power(domain.Count, b) / power(domain.Count - a, b - a_b);
        p2 := log10(p2);
        Domain_weighted := p1 + p2;

      end;
      if max_dom < Domain_weighted then
        max_dom := Domain_weighted;
      if min_dom > Domain_weighted then
        min_dom := Domain_weighted;
      Matrix_Edge[i][j][1] := Domain_weighted;
      Matrix_Edge[j][i][1] := Matrix_Edge[i][j][1];

    //----------------------complex---------------------------------------------
      if (NodeInfo[i].ComplexNum * NodeInfo[j].ComplexNum > 0) then
      begin
        Complex_Weight := 0;
        a := Node_Info[i].ComplexNum;
        b := Node_Info[j].ComplexNum;
        a_b := 0;
        for k := 0 to high(Detect_Node) do
          a_b := a_b + Matrix_Complex[Neighbor_list[i][l]][k] * Matrix_Complex[Neighbor_list[j][m]][k];

        Complex_Weight := a_b / sqrt(a * b);
        Matrix_Edge[i][j][2] := Complex_Weight;
        Matrix_Edge[j][i][2] := Complex_Weight;
      end;
    end;
    application.ProcessMessages;
  end;
  for i := 0 to nodes.Count - 1 do
    for j := i + 1 to nodes.Count - 1 do
    begin
      Matrix_Edge[i][j][1] := (Matrix_Edge[i][j][1] - min_dom) / (max_dom - min_dom);
      Matrix_Edge[j][i][1] := Matrix_Edge[i][j][1];
    end;
  //setlength(MatrixEdge, 0);
  setlength(Edge_Info, 0);
  setlength(NodeInfo, 0);
  //--------------------------------
  setlength(Valid_Proteins, 0);
  for i := low(Node_Info) to high(Node_Info) do
  begin
    if Node_Info[i].GoNums = 0 then
      continue;
    setlength(Valid_Proteins, length(Valid_Proteins) + 1);
    Valid_Proteins[high(Valid_Proteins)] := i;
  end;

  if CheckBox1.Checked = true then
  begin
    L := Round(length(Node_Info) * StrToFloat(Edit1.Text));
    Setlength(Select_Proteins, 0);
    for i := 0 to L - 1 do
    begin
      Weighted := False;
      while Weighted = False do
      begin
        Randomize;
        IPos := Random(length(Node_Info));
        for j := 0 to high(Select_Proteins) do
          if Select_Proteins[j] = IPos then
            break;
        if j > high(Select_Proteins) then
          Weighted := True;
      end;
      setlength(Select_Proteins, length(Select_Proteins) + 1);
      Select_Proteins[high(Select_Proteins)] := IPos;
     // Node_Info[IPos].Unknown := 1;
      Application.ProcessMessages;
    end;
  end;
  //----------------------------------------------------------------------------
  Btn_Init.Enabled := True;
  Btn_NGF.Enabled := True;
end;
2. Construction of the propagation network
procedure TFrmMain.Btn_DistanceClick(Sender: TObject);
var
  i, j, k, l, h, In_Num, Out_Num, temp2, MatchNum, A_Index, Common, Common1, A1, B1: integer;
  New_pr, New_newpr, h0, vd, vc, vd0, vc0: array of single;
  a, standard, Temp, FuncProtein, MaxSim, Precision, Value: double;
  Matrix_A: array of array of single;
  firsttime, lasttime: real;
begin
  Btn_Distance.Enabled := false;
  Gauge1.MaxValue := length(Valid_Proteins);
  Gauge1.MinValue := 0;
  Gauge1.Progress := 0;
  //----------------------prepare-----------------------------
  a := strtofloat(edit1.Text);
  setlength(New_pr, length(Valid_Proteins));
  setlength(New_newpr, length(Valid_Proteins));
  setlength(vd, length(Valid_Proteins));
  setlength(vd0, length(Valid_Proteins));
  setlength(vc, length(Valid_Proteins));
  setlength(vc0, length(Valid_Proteins));
  setlength(New_Matrix, Nodes.Count);
  setlength(Matrix_A, length(Valid_Proteins));
  //----------------------Initial----------------------------
  for i := 0 to nodes.Count - 1 do
    setlength(New_Matrix[i], Nodes.Count);
  for i := low(Valid_Proteins) to high(Valid_Proteins) do
  begin
    setlength(Matrix_A[i], length(Valid_Proteins));
    temp := 0;
    for j := low(Valid_Proteins) to high(Valid_Proteins) do
      temp := temp + Matrix_Edge[Valid_Proteins[i]][Valid_Proteins[j]][0];
    for j := low(Valid_Proteins) to high(Valid_Proteins) do
      if temp <> 0 then
        Matrix_A[i][j] := (Matrix_Edge[Valid_Proteins[i]][Valid_Proteins[j]][0]) / temp
      else
        Matrix_A[i][j] := 0;
  end;
 //---------------------Random walk-----------------------------------
  for i := low(Valid_Proteins) to high(Valid_Proteins) do //random walk for each node
  begin
    Gauge1.Progress := Gauge1.Progress + 1;
    for j := low(Valid_Proteins) to high(Valid_Proteins) do
    begin
      vd0[j] := Matrix_Edge[Valid_Proteins[i]][Valid_Proteins[j]][1];
      vc0[j] := Matrix_Edge[Valid_Proteins[i]][Valid_Proteins[j]][2];
      vc[j] := 1.0 / length(Valid_Proteins);
      vd[j] := 1.0 / length(Valid_Proteins);
      New_pr[j] := 0;
    end;
    for k := 0 to 19 do //iteration
    begin
      for l := low(Valid_Proteins) to high(Valid_Proteins) do
      begin
        Temp := 0;
        for j := low(Valid_Proteins) to high(Valid_Proteins) do
          Temp := Temp + Matrix_A[l][j] * vc[j];
        vd[l] := (1 - a) * vd0[l] + a * Temp;
        Temp := 0;
        for j := low(Valid_Proteins) to high(Valid_Proteins) do
          Temp := Temp + Matrix_A[l][j] * vd[j];
        vc[l] := (1 - a) * vc0[l] + a * Temp;
        New_newpr[l] := (vd[l] + vc[l]);
        application.ProcessMessages;
      end;
      Temp := 0;
      for l := low(Valid_Proteins) to high(Valid_Proteins) do
        Temp := Temp + abs(New_pr[l] - New_newpr[l]);
      if temp < 0.001 then
        break;

      for l := low(Valid_Proteins) to high(Valid_Proteins) do
        New_pr[l] := New_newpr[l];
    end; //for k := 0 to 19 do //iteration
    //Memo4.Lines.Add(IntTostr(k));
    for j := low(Valid_Proteins) to high(Valid_Proteins) do
      if New_newpr[j] > 0.001 then
        New_Matrix[Valid_Proteins[i]][Valid_Proteins[j]] := New_newpr[j]
      else
        New_Matrix[Valid_Proteins[i]][Valid_Proteins[j]] := 0;
  end; //for i := low(Valid_Proteins) to high(Valid_Proteins) do
  for i := 0 to nodes.Count - 1 do
    for j := i + 1 to nodes.Count - 1 do
      if New_Matrix[i][j] * New_Matrix[j][i] = 0 then
      begin
        New_Matrix[i][j] := 0;
        New_Matrix[j][i] := 0;
      end
      else
      begin
        New_Matrix[i][j] := (New_Matrix[i][j] + New_Matrix[j][i]) / 2;
        New_Matrix[i][j] := New_Matrix[j][i];
      end;
  Btn_Distance.Enabled := true;
end;

3. NGF algorithm
procedure TFrmMain.Btn_NGFClick(Sender: TObject);
var
  i, j, k, l, m, n, p, selected, temp, Certain_Go_Num, Nums, MaxRange, ipos, TotalNum, BenchBark_functions, Predict_Functions, Protein_Num, Predicted_Num: integer;
  candidate_protein, Candidate_GO: array of integer;
  candidate_protein_score, Go_Score, GO_Nums: array of double;
  Value, Precision, Recall, TotalValue, MaxVlue, Current_Density, Temp_Density, In_Density, out_Density: Double;
  Candidate_Index, Neighbor_Set, Temp_set: array of integer; //候选节点集合
  pre, rec: array[0..0] of double;
  TP, FP: array[0..0] of integer;
  firsttime, lasttime: real;
  predictedFunctions, s: string;
  Existed: boolean;
begin
  firsttime := GetTickCount;
  Btn_NGF.Enabled := False;
  Btn_Distance.Click;
  Precision := 0; Recall := 0;
  for i := 0 to 0 do
  begin
    pre[i] := 0;
    rec[i] := 0;
    FP[i] := 0;
    TP[i] := 0;
  end;
  TotalNum := 0;
  Protein_Num := 0;
  BenchBark_functions := 0;
  Predict_Functions := 0;
  Gauge1.MaxValue := length(Valid_Proteins);
  Gauge1.MinValue := 0;
  Gauge1.Progress := 0;
           //-----------------------------------------------------------------------------
  setlength(Protein_list, length(Valid_Proteins));
  for i := low(Valid_Proteins) to high(Valid_Proteins) do
    for j := i + 1 to high(Valid_Proteins) do
      if New_Matrix[Valid_Proteins[i]][Valid_Proteins[j]] <> 0 then
      begin
        setlength(Protein_list[i], length(Protein_list[i]) + 1);
        Protein_list[i][high(Protein_list[i])] := Valid_Proteins[j];
        setlength(Protein_list[j], length(Protein_list[j]) + 1);
        Protein_list[j][high(Protein_list[j])] := Valid_Proteins[i];
      end;
  //--------------------------------------------------------------------------
  memo3.Lines.Clear;
  Memo4.Lines.Clear;
  Predicted_Num := 0;
  for i := low(Valid_Proteins) to high(Valid_Proteins) do
  begin
    if (CheckBox1.Checked = True) and (Node_Info[Valid_Proteins[i]].Unknown = 0) then continue;
    Gauge1.Progress := Gauge1.Progress + 1;
    Inc(Predicted_Num);
    s := Node_Info[Valid_Proteins[i]].Node;
    Inc(Protein_Num);
    predictedFunctions := '';
    BenchBark_functions := BenchBark_functions + Node_Info[Valid_Proteins[i]].GONums;
   //------------------------------------------------------------
    setlength(Candidate_GO, 0);
    setlength(GO_Score, 0);
    setlength(GO_Nums, 0);
    MaxRange := 0;
    MaxVlue := 0; ipos := -1;
    setlength(candidate_protein, 0);
    setlength(candidate_protein_score, 0);
    m := 1;
      //------------------------------------------------
    for j := 0 to high(protein_list[i]) do
    begin
      if Node_Info[protein_list[i][j]].GoNums = 0 then continue;
      if (CheckBox1.Checked = True) and (Node_Info[protein_list[i][j]].Unknown = 1) then continue;
      if New_Matrix[Valid_Proteins[i]][protein_list[i][j]] <> 0 then
      begin
        In_Queue(protein_list[i][j], New_Matrix[Valid_Proteins[i]][protein_list[i][j]]);
        m := m + 1;
      end;
    end;
                           //--------------------------高内聚-----------------------------------------
    Selected := Out_Queue(Value);
    setlength(Candidate_Index, 1);
    setlength(Neighbor_Set, 0);
    Candidate_Index[0] := Valid_Proteins[i];
    Current_Density := 0;
    while Selected >= 0 do
    begin
      //-------------------计算节点加入之前的适用度-----------------------------
      In_Density := 0;
      for j := 0 to high(Candidate_Index) do
        for k := j + 1 to high(Candidate_Index) do
          In_Density := In_Density + New_Matrix[Candidate_Index[j]][Candidate_Index[k]];
      In_Density := In_Density * 2;
      out_Density := 0;
      for j := 0 to high(Valid_Proteins) do
      begin
        for k := 0 to high(Candidate_Index) do
          if (Candidate_Index[k] = Valid_Proteins[j]) then
            break;
        if k <= high(Candidate_Index) then
          continue;
        for k := 0 to high(Candidate_Index) do
          out_Density := out_Density + New_Matrix[Candidate_Index[k]][Valid_Proteins[j]];
      end;
      out_Density := out_Density * 2;
      Temp_Density := In_Density / (In_Density + out_Density);
      //------------------------------------------------
      setlength(Candidate_Index, length(Candidate_Index) + 1);
      Candidate_Index[high(Candidate_Index)] := Selected;
      In_Density := 0;
      for j := 0 to high(Candidate_Index) do
        for k := j + 1 to high(Candidate_Index) do
          In_Density := In_Density + New_Matrix[Candidate_Index[j]][Candidate_Index[k]];
      In_Density := In_Density * 2;
      out_Density := 0;
      for j := 0 to high(Valid_Proteins) do
      begin
        for k := 0 to high(Candidate_Index) do
          if (Candidate_Index[k] = Valid_Proteins[j]) then
            break;
        if k <= high(Candidate_Index) then
          continue;
        for k := 0 to high(Candidate_Index) do
          out_Density := out_Density + New_Matrix[Candidate_Index[k]][Valid_Proteins[j]];
      end;
      out_Density := out_Density * 2;
      Current_Density := In_Density / (In_Density + out_Density);
      if Current_Density - Temp_Density < 0 then
        setlength(Candidate_Index, length(Candidate_Index) - 1); //移除这个节点
      Selected := Out_Queue(Value);
      Application.ProcessMessages;
    end; //while Selected >= 0 do
           //-------------------------------------------------------------------------
    if length(Candidate_Index) > 1 then
    begin
      setlength(candidate_protein, length(Candidate_Index) - 1);
      setlength(candidate_protein_score, length(Candidate_Index) - 1);
      s := '';
      for k := 1 to high(Candidate_Index) do
      begin
        s := s + Node_Info[Candidate_Index[k]].Node;
        if k < high(Candidate_Index) then
          s := s + ',';
        candidate_protein[k - 1] := Candidate_Index[k];
        candidate_protein_score[k - 1] := New_Matrix[Valid_Proteins[i]][Candidate_Index[k]];
        if MaxVlue < New_Matrix[Valid_Proteins[i]][Candidate_Index[k]] then
        begin
          MaxVlue := New_Matrix[Valid_Proteins[i]][Candidate_Index[k]];
          ipos := Candidate_Index[k];
        end;
      end;
      Memo4.Lines.Add(Node_INfo[Valid_Proteins[i]].Node + #9 + IntToStr(Node_INfo[Valid_Proteins[i]].GoNums) + #9 + IntToStr(length(Candidate_Index)) + #9 + FloatToStr(Current_Density) + #9 + s);
    end;
                    //-------------------------------------------------------------------------
    for j := 0 to Go.Count - 1 do
    begin
      Value := 0; Certain_Go_Num := 0;
      for n := 0 to length(candidate_protein) - 1 do
        if Matrix_GO[candidate_protein[n]][j] = 1 then
        begin
          Value := Value + candidate_protein_score[n];
          Inc(Certain_Go_Num);
        end;
      if Certain_Go_Num = 0 then
        continue;
      setlength(Candidate_GO, length(Candidate_GO) + 1);
      Candidate_GO[high(Candidate_GO)] := j;
      setlength(GO_Score, length(GO_Score) + 1);
      GO_Score[high(GO_Score)] := Value;
      setlength(GO_Nums, length(GO_Nums) + 1);
      GO_Nums[high(GO_Nums)] := 1;
    end;
    if MaxRange = 0 then
      MaxRange := Node_Info[ipos].GONums;
   //---------------------- --------------------------------------
    for j := low(Candidate_GO) to high(Candidate_GO) - 1 do
      for k := j + 1 to high(Candidate_GO) do
        if GO_Score[j] < GO_Score[k] then
        begin
          Value := GO_Score[j];
          GO_Score[j] := GO_Score[k];
          GO_Score[k] := Value;
          m := Candidate_GO[j];
          Candidate_GO[j] := Candidate_GO[k];
          Candidate_GO[k] := m;
        end;
   //--------------------- ----------------------------------------------
   //MaxRange := Node_Info[Valid_Proteins[i]].GONums;
    for l := 0 to high(Pre) do
    begin
      if length(Pre) > 1 then
      begin
        MaxRange := length(Candidate_GO);
        if MaxRange > l + 1 then
          MaxRange := l + 1;
      end
      else if length(Candidate_GO) < MaxRange then
        MaxRange := length(Candidate_GO);

      Nums := 0;
      for k := 0 to MaxRange - 1 do
      begin
        for j := 0 to Go.Count - 1 do
        begin
          if Matrix_GO[Valid_Proteins[i]][j] = 0 then
            continue;
          if Candidate_GO[k] = j then
          begin
            Inc(Nums);
            break;
          end;
        end;
       
 predictedFunctions := predictedFunctions + GO.Strings[Candidate_GO[k]];
        if j = GO.Count then
          predictedFunctions := predictedFunctions + '(F)';
        predictedFunctions := predictedFunctions + ' ';
      end;

      if Nums <> 0 then
      begin
        Inc(TotalNum);
        Predict_Functions := Predict_Functions + Nums;
      end;
      if MaxRange <> 0 then
        Pre[l] := Pre[l] + 1.0 * Nums / MaxRange;
      Rec[l] := Rec[l] + 1.0 * Nums / (Node_Info[Valid_Proteins[i]].GONums);
      TP[l] := TP[l] + Nums;
      FP[l] := FP[l] + (MaxRange - Nums);
    end; //for l := 0 to high(Prec) do

    Memo3.Lines.Add(Node_Info[Valid_Proteins[i]].Node + #9 + IntToStr(Nums) + #9 + FloatToStr(MaxRange) + #9 + IntToStr(Node_Info[Valid_Proteins[i]].GONums));
   //---------------------------------------------------------------------------
    Application.ProcessMessages;
  end; //for i := low(Valid_Proteins) to high(Valid_Proteins) do
 //-----------------------------------------------------------------------------
  Btn_NGF.Enabled := True;


end;
