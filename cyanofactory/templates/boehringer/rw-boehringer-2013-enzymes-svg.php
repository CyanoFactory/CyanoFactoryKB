<html>
<head>
<title>Boehringer Pathways</title>

    <style type="text/css">
        #scrolly{
            overflow-x: auto;
            overflow-y: auto;
            width: 85%;
            height: 85%;
            margin: auto;
            border: 5px;
            border: 1px solid;
            border-radius: 25px;
            box-shadow: 10px 10px 5px #888888;
            white-space: nowrap;
            }
        #sidebar{
        font-family:"Verdana", Sans-serif;
        font-size: 13px;
        width: 100px;
        height: 85%;
        margin: 50px 20px;
        border: 5px;
        border: 1px solid;
        border-radius: 15px;
        box-shadow: 10px 10px 5px #888888;
        float: right
        }
        h2{
        text-align:center;
        color:blue;
        font-family:"Verdana", Sans-serif
        }
        h5{
        text-align:center;
        color:black;
        font-family:"Verdana", Sans-serif
        }
		.roundrect{
		border-radius: 10px;
		}
    </style>

</head>
<body>


<div id="sidebar"  align="left">
<img class="roundrect" src="logo-facebook.png" width="100"/>
<br>
<br>
<form action="rw-boehringer-2013-enzymes-svg.php" method="post" >
  <br><br><p style="margin-left: 3" align="left"><b>ECs:</b>  Paste white-space separated EC numbers and/or keywords here. Both can be linked to colors with #.<br><br>
  <?php
  if($_REQUEST["request"] == "") {
  	$req_form = "1.1.1.1 2.2.2.2 4.1.2.20#green 1.1.1.20 1.2.1.3 ascorbate#red";
  	}
  	else{
  	$req_form = $_REQUEST["request"];
  	
  	}
  ?>
  <textarea name="request" cols="10" rows="10"><?php echo $req_form; ?></textarea>
  <br>
  <input type="submit">
  </p>
</form>

</div>


<h2> CyanoFactory - Annotated Pathway Map </h2>

<!-- viewbox -->

<div id="scrolly" >

<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="6871" height="4804" pointer-events="visible"> 

<image x="0" y="0" width="6871" height="4804" xlink:href="boepathways.gif" />

<?php

//print $_POST["request"];
// DATABASE
// echo sqlite_libversion();
// echo phpversion();

$db = new SQLite3('pathway.db');

$lines = preg_split("/[\s\n\r]/", $_REQUEST["request"]); // ?request=4.1.2.20-orange%201.2.1.3-green

foreach ($lines as $i){
//	echo "<b>".$i."</b><br>";
	$line = preg_split("/#/", $i);
	//echo "<!--".$line[0]."-->\n";
	if($line[0] == ""){continue;}
	
	if(preg_match("/[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+/",$line[0])){
	// IS EC
		$query = "SELECT * FROM enzymes WHERE ec=\"".$line[0]."\"";
		//echo "<!-- ".$query." -->\n";
		$results = $db->query($query);
		//if($results == FALSE){echo "<!-- NOO HIT -->\n";}
		while ($row = $results->fetchArray(SQLITE3_ASSOC)) {
   		 	//var_dump($row);
    		//echo $row["ec"]."-".$row["color"]."<br>";
    		//if(!isset($row["ec"])){echo "<!-- NOO HIT -->\n";};
	    	$eccount++;
    		//echo "<!-- EC ".$row["ec"]." -->\n";
  			if($line[1] != ""){$row["color"] = $line[1];}
			echo "<a xlink:href=\"http://www.brenda-enzymes.org/php/result_flat.php4?ecno=".$row["ec"]."\" xlink:title=\"".$row["title"]."\">";
			echo "\n";
			echo '<rect x="'.$row["x"].'" y="'.$row["y"].'" width="'.$row["w"].'" height="'.$row["h"].'" style="fill:'.$row["color"].'; fill-opacity:0.2"/>';
			echo "\n";
			echo "</a>";
			echo "\n";
			}
		if ($eccount==0){echo "<!-- NO HIT for ".$line[0]." -->\n"; $eccount=0; $ecnohit++;}else{$eccount=0; $echit++;}
	
	}
	
	else{
		// IS NO EC
		$query = "SELECT * FROM metabolites WHERE title LIKE \"%".$line[0]."%\"";
		$results = $db->query($query);
		while ($row = $results->fetchArray(SQLITE3_ASSOC)) {
			$metcount++;
			if($line[1] != ""){$row["color"] = $line[1];}
			echo '<rect x="'.$row["x"].'" y="'.$row["y"].'" width="'.$row["w"].'" height="'.$row["h"].'" style="fill:'.$row["color"].'; fill-opacity:0.2"/>';
			echo "\n";
		}
		if ($metcount==0){echo "<!-- NO HIT for ".$line[0]." -->\n"; $metcount=0; $metnohit++;}else{$metcount=0; $methit++;}
	}
	
	
	}


/* METABOLITES
$query = "SELECT * FROM metabolites";
$results = $db->query($query);
while ($row = $results->fetchArray(SQLITE3_ASSOC)) {
	echo '<rect x="'.$row["x"].'" y="'.$row["y"].'" width="'.$row["w"].'" height="'.$row["h"].'" style="fill:'.$row["color"].'; fill-opacity:0.2"/>';
	echo "\n";
	}
*/
?>
   

</svg>
<!-- <img src="boepathways.gif" border="0" height="4804" wigth="6871"/> -->
</div>

<h5> by Robbe Wunschiers / Version 2013.5 / ECs: <?php echo "hits=".$echit." - missing=".$ecnohit; ?> / Metabolites: <?php echo "hits=".$methit." - missing=".$metnohit; ?></h5>

</body>
</html>
