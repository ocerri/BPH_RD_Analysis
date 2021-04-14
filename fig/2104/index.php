<html>
<head>
<title><?php echo getcwd(); ?></title>
<style type='text/css'>
body {
    font-family: "Candara", sans-serif;
    font-size: 9pt;
    line-height: 10.5pt;
}
div.pic h3 {
    font-size: 11pt;
    margin: 0.5em 1em 0.2em 1em;
}
div.pic p {
    font-size: 11pt;
    margin: 0.2em 1em 0.1em 1em;
}
div.pic {
    display: block;
    float: left;
    background-color: white;
    border: 1px solid #ccc;
    padding: 2px;
    text-align: center;
    margin: 5px 1%;
    -moz-box-shadow: 7px 5px 5px rgb(80,80,80);    /* Firefox 3.5 */
    -webkit-box-shadow: 7px 5px 5px rgb(80,80,80); /* Chrome, Safari */
    box-shadow: 7px 5px 5px rgb(80,80,80);         /* New browsers */
}
a { text-decoration: none; color: rgb(80,0,0); }
a:hover { text-decoration: underline; color: rgb(255,80,80); }
</style>
</head>
<body>
<?php
$match = htmlspecialchars($_GET['match']);
$noplots = htmlspecialchars($_GET['noplots']);
?>
<h1><?php echo getcwd(); ?></h1>
<div style="display: block; clear:both;">
<p><form>Filter: <input type="text" name="match" size="30" value="<?php if (isset($match)) print $match;  ?>" /><input type="Submit" value="Go" /></form></p>
<h2><a name="directories">Directories</a></h2>
<ul>
<?php
foreach (glob("*") as $filename) {
    if ($noplots || !in_array($filename, $displayed)) {
        if (isset($match) && !fnmatch('*'.$match.'*', $filename)) continue;
        if (is_dir($filename)) {
          print "<li>[DIR] <a href=\"$filename\">$filename</a></li>";
        }
    }
}
?>
</ul>
</div>
<h2><a name="plots">Plots</a></h2>
<div style="width: 100%;">
<?php
$displayed = array();
if ($noplots) {
    print "Plots will not be displayed.\n";
} else {
    $other_exts = array('.pdf', '.cxx', '.eps', '.root', '.txt');
    $filenames = glob("*.png"); sort($filenames);
    foreach ($filenames as $filename) {
        if (isset($match) && !fnmatch('*'.$match.'*', $filename)) continue;
        array_push($displayed, $filename);
        print "<div class='pic' style=\"width: 30%;\">\n";
        print "<h3><a href=\"$filename\">$filename</a></h3>";
        print "<a href=\"$filename\"><img src=\"$filename\" style=\"border: none; width: 100%; \"></a>";
        $others = array();
        foreach ($other_exts as $ex) {
            $other_filename = str_replace('.png', $ex, $filename);
            if (file_exists($other_filename)) {
                array_push($others, "<a class=\"file\" href=\"$other_filename\">[" . $ex . "]</a>");
                if ($ex != '.txt') array_push($displayed, $other_filename);
            }
        }
        if ($others) print "<p>Also as ".implode(', ',$others)."</p>";
        print "</div>";
    }
}
?>
</div>
<div style="display: block; clear:both;">
<h2><a name="files">Other files</a></h2>
<ul>
<?php
foreach (glob("*") as $filename) {
    if ($noplots || !in_array($filename, $displayed)) {
        if (isset($match) && !fnmatch('*'.$match.'*', $filename)) continue;
        if (!is_dir($filename)) {
            print "<li><a href=\"$filename\">$filename</a></li>";
        }
    }
}
?>
</ul>
</div>
</body>
</html>
