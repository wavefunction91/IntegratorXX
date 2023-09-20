#pragma once

namespace IntegratorXX {
namespace AhrensBeylkinGrids {

/**
 *  \brief Ahrens-Beylkin Quadrature specification for order 29 grid with 312 points
 * 
 */
template <typename T>
struct ahrens_beylkin_312 {

  static constexpr std::array<cartesian_pt_t<T>,312> points = {
      0.000000000000000E+00,         -0.5257311121191336,          0.8506508083520399,
      0.000000000000000E+00,         -0.5257311121191336,         -0.8506508083520399,
         0.8506508083520399,       0.000000000000000E+00,         -0.5257311121191336,
         0.8506508083520399,       0.000000000000000E+00,          0.5257311121191336,
        -0.5257311121191336,          0.8506508083520399,       0.000000000000000E+00,
        -0.5257311121191336,         -0.8506508083520399,       0.000000000000000E+00,
      0.000000000000000E+00,          0.5257311121191336,         -0.8506508083520399,
      0.000000000000000E+00,          0.5257311121191336,          0.8506508083520399,
        -0.8506508083520399,       0.000000000000000E+00,          0.5257311121191336,
        -0.8506508083520399,       0.000000000000000E+00,         -0.5257311121191336,
         0.5257311121191336,         -0.8506508083520399,       0.000000000000000E+00,
         0.5257311121191336,          0.8506508083520399,       0.000000000000000E+00,
        -0.1201736413467416,         -0.9416989679914392,          0.3142631852593033,
        -0.9416989679914392,         -0.3142631852593033,          0.1201736413467416,
         0.5675831786087092,          0.6651846670726805,         -0.4851584216025124,
        -0.9561177587721386,         -0.2765143009187586,      -9.684354778689629E-02,
        -0.1201736413467416,          0.9416989679914392,         -0.3142631852593033,
         0.9416989679914392,         -0.3142631852593033,         -0.1201736413467416,
         0.6053320629492540,         -0.6418545735128353,         -0.4707396308218130,
     -9.684354778689629E-02,         -0.9561177587721386,         -0.2765143009187586,
         0.2765143009187586,      -9.684354778689629E-02,          0.9561177587721386,
        -0.4707396308218130,         -0.6053320629492540,          0.6418545735128353,
        -0.6651846670726805,          0.4851584216025124,          0.5675831786087092,
         0.4709593371696262,      -2.333009355984531E-02,          0.8818463638680126,
         0.4851584216025124,          0.5675831786087092,         -0.6651846670726805,
      2.333009355984531E-02,          0.8818463638680126,         -0.4709593371696262,
         0.6418545735128353,          0.4707396308218130,          0.6053320629492540,
        -0.8818463638680126,         -0.4709593371696262,      -2.333009355984531E-02,
        -0.9561177587721386,          0.2765143009187586,       9.684354778689629E-02,
        -0.3142631852593033,          0.1201736413467416,         -0.9416989679914392,
        -0.3142631852593033,         -0.1201736413467416,          0.9416989679914392,
         0.5675831786087092,         -0.6651846670726805,          0.4851584216025124,
         0.3142631852593033,         -0.1201736413467416,         -0.9416989679914392,
        -0.2765143009187586,      -9.684354778689629E-02,         -0.9561177587721386,
         0.4707396308218130,         -0.6053320629492540,         -0.6418545735128353,
        -0.4709593371696262,      -2.333009355984531E-02,         -0.8818463638680126,
      9.684354778689629E-02,          0.9561177587721386,         -0.2765143009187586,
        -0.4851584216025124,         -0.5675831786087092,         -0.6651846670726805,
        -0.6053320629492540,          0.6418545735128353,         -0.4707396308218130,
     -2.333009355984531E-02,         -0.8818463638680126,         -0.4709593371696262,
         0.6651846670726805,          0.4851584216025124,         -0.5675831786087092,
         0.1201736413467416,          0.9416989679914392,          0.3142631852593033,
        -0.9416989679914392,          0.3142631852593033,         -0.1201736413467416,
        -0.8818463638680126,          0.4709593371696262,       2.333009355984531E-02,
         0.6418545735128353,         -0.4707396308218130,         -0.6053320629492540,
         0.1201736413467416,         -0.9416989679914392,         -0.3142631852593033,
         0.9416989679914392,          0.3142631852593033,          0.1201736413467416,
         0.3142631852593033,          0.1201736413467416,          0.9416989679914392,
         0.4707396308218130,          0.6053320629492540,          0.6418545735128353,
         0.6651846670726805,         -0.4851584216025124,          0.5675831786087092,
        -0.4709593371696262,       2.333009355984531E-02,          0.8818463638680126,
        -0.2765143009187586,       9.684354778689629E-02,          0.9561177587721386,
     -9.684354778689629E-02,          0.9561177587721386,          0.2765143009187586,
         0.6053320629492540,          0.6418545735128353,          0.4707396308218130,
         0.9561177587721386,         -0.2765143009187586,       9.684354778689629E-02,
        -0.6418545735128353,         -0.4707396308218130,          0.6053320629492540,
         0.8818463638680126,          0.4709593371696262,      -2.333009355984531E-02,
        -0.5675831786087092,          0.6651846670726805,          0.4851584216025124,
         0.4709593371696262,       2.333009355984531E-02,         -0.8818463638680126,
        -0.6651846670726805,         -0.4851584216025124,         -0.5675831786087092,
     -2.333009355984531E-02,          0.8818463638680126,          0.4709593371696262,
      9.684354778689629E-02,         -0.9561177587721386,          0.2765143009187586,
        -0.6418545735128353,          0.4707396308218130,         -0.6053320629492540,
         0.9561177587721386,          0.2765143009187586,      -9.684354778689629E-02,
        -0.4707396308218130,          0.6053320629492540,         -0.6418545735128353,
         0.2765143009187586,       9.684354778689629E-02,         -0.9561177587721386,
        -0.4851584216025124,          0.5675831786087092,          0.6651846670726805,
        -0.6053320629492540,         -0.6418545735128353,          0.4707396308218130,
        -0.5675831786087092,         -0.6651846670726805,         -0.4851584216025124,
         0.8818463638680126,         -0.4709593371696262,       2.333009355984531E-02,
      2.333009355984531E-02,         -0.8818463638680126,          0.4709593371696262,
         0.4851584216025124,         -0.5675831786087092,          0.6651846670726805,
        -0.4470439773240164,         -0.8796732818399999,          0.1622547366185910,
        -0.8796732818399999,         -0.1622547366185910,          0.4470439773240164,
         0.4923990799707575,          0.8516422868410736,         -0.1795788442664059,
        -0.9309421890415280,      -2.803099499892637E-02,         -0.3640891429058797,
        -0.4470439773240164,          0.8796732818399999,         -0.1622547366185910,
         0.8796732818399999,         -0.1622547366185910,         -0.4470439773240164,
         0.6266228215904223,         -0.7686874524229369,         -0.1283099370648778,
        -0.3640891429058797,         -0.9309421890415280,      -2.803099499892637E-02,
      2.803099499892637E-02,         -0.3640891429058797,          0.9309421890415280,
        -0.1283099370648778,         -0.6266228215904223,          0.7686874524229369,
        -0.8516422868410736,          0.1795788442664059,          0.4923990799707575,
         0.7513633447751221,      -8.295483441813664E-02,          0.6546538165893486,
         0.1795788442664059,          0.4923990799707575,         -0.8516422868410736,
      8.295483441813664E-02,          0.6546538165893486,         -0.7513633447751221,
         0.7686874524229369,          0.1283099370648778,          0.6266228215904223,
        -0.6546538165893486,         -0.7513633447751221,      -8.295483441813664E-02,
        -0.9309421890415280,       2.803099499892637E-02,          0.3640891429058797,
        -0.1622547366185910,          0.4470439773240164,         -0.8796732818399999,
        -0.1622547366185910,         -0.4470439773240164,          0.8796732818399999,
         0.4923990799707575,         -0.8516422868410736,          0.1795788442664059,
         0.1622547366185910,         -0.4470439773240164,         -0.8796732818399999,
     -2.803099499892637E-02,         -0.3640891429058797,         -0.9309421890415280,
         0.1283099370648778,         -0.6266228215904223,         -0.7686874524229369,
        -0.7513633447751221,      -8.295483441813664E-02,         -0.6546538165893486,
         0.3640891429058797,          0.9309421890415280,      -2.803099499892637E-02,
        -0.1795788442664059,         -0.4923990799707575,         -0.8516422868410736,
        -0.6266228215904223,          0.7686874524229369,         -0.1283099370648778,
     -8.295483441813664E-02,         -0.6546538165893486,         -0.7513633447751221,
         0.8516422868410736,          0.1795788442664059,         -0.4923990799707575,
         0.4470439773240164,          0.8796732818399999,          0.1622547366185910,
        -0.8796732818399999,          0.1622547366185910,         -0.4470439773240164,
        -0.6546538165893486,          0.7513633447751221,       8.295483441813664E-02,
         0.7686874524229369,         -0.1283099370648778,         -0.6266228215904223,
         0.4470439773240164,         -0.8796732818399999,         -0.1622547366185910,
         0.8796732818399999,          0.1622547366185910,          0.4470439773240164,
         0.1622547366185910,          0.4470439773240164,          0.8796732818399999,
         0.1283099370648778,          0.6266228215904223,          0.7686874524229369,
         0.8516422868410736,         -0.1795788442664059,          0.4923990799707575,
        -0.7513633447751221,       8.295483441813664E-02,          0.6546538165893486,
     -2.803099499892637E-02,          0.3640891429058797,          0.9309421890415280,
        -0.3640891429058797,          0.9309421890415280,       2.803099499892637E-02,
         0.6266228215904223,          0.7686874524229369,          0.1283099370648778,
         0.9309421890415280,      -2.803099499892637E-02,          0.3640891429058797,
        -0.7686874524229369,         -0.1283099370648778,          0.6266228215904223,
         0.6546538165893486,          0.7513633447751221,      -8.295483441813664E-02,
        -0.4923990799707575,          0.8516422868410736,          0.1795788442664059,
         0.7513633447751221,       8.295483441813664E-02,         -0.6546538165893486,
        -0.8516422868410736,         -0.1795788442664059,         -0.4923990799707575,
     -8.295483441813664E-02,          0.6546538165893486,          0.7513633447751221,
         0.3640891429058797,         -0.9309421890415280,       2.803099499892637E-02,
        -0.7686874524229369,          0.1283099370648778,         -0.6266228215904223,
         0.9309421890415280,       2.803099499892637E-02,         -0.3640891429058797,
        -0.1283099370648778,          0.6266228215904223,         -0.7686874524229369,
      2.803099499892637E-02,          0.3640891429058797,         -0.9309421890415280,
        -0.1795788442664059,          0.4923990799707575,          0.8516422868410736,
        -0.6266228215904223,         -0.7686874524229369,          0.1283099370648778,
        -0.4923990799707575,         -0.8516422868410736,         -0.1795788442664059,
         0.6546538165893486,         -0.7513633447751221,       8.295483441813664E-02,
      8.295483441813664E-02,         -0.6546538165893486,          0.7513633447751221,
         0.1795788442664059,         -0.4923990799707575,          0.8516422868410736,
        -0.3358108236422033,         -0.8807282002850116,          0.3339894129272207,
        -0.8807282002850116,         -0.3339894129272207,          0.3358108236422033,
         0.4417581235917372,          0.8152491678999195,         -0.3744576804916627,
        -0.9832900393199363,      -6.547903238509202E-02,         -0.1698622821349993,
        -0.3358108236422033,          0.8807282002850116,         -0.3339894129272207,
         0.8807282002850116,         -0.3339894129272207,         -0.3358108236422033,
         0.7102685041338660,         -0.6493006263927156,         -0.2718958414567380,
        -0.1698622821349993,         -0.9832900393199363,      -6.547903238509202E-02,
      6.547903238509202E-02,         -0.1698622821349993,          0.9832900393199363,
        -0.2718958414567380,         -0.7102685041338660,          0.6493006263927156,
        -0.8152491678999195,          0.3744576804916627,          0.4417581235917372,
         0.6088323588282737,         -0.1659485415072040,          0.7757475365189580,
         0.3744576804916627,          0.4417581235917372,         -0.8152491678999195,
         0.1659485415072040,          0.7757475365189580,         -0.6088323588282737,
         0.6493006263927156,          0.2718958414567380,          0.7102685041338660,
        -0.7757475365189580,         -0.6088323588282737,         -0.1659485415072040,
        -0.9832900393199363,       6.547903238509202E-02,          0.1698622821349993,
        -0.3339894129272207,          0.3358108236422033,         -0.8807282002850116,
        -0.3339894129272207,         -0.3358108236422033,          0.8807282002850116,
         0.4417581235917372,         -0.8152491678999195,          0.3744576804916627,
         0.3339894129272207,         -0.3358108236422033,         -0.8807282002850116,
     -6.547903238509202E-02,         -0.1698622821349993,         -0.9832900393199363,
         0.2718958414567380,         -0.7102685041338660,         -0.6493006263927156,
        -0.6088323588282737,         -0.1659485415072040,         -0.7757475365189580,
         0.1698622821349993,          0.9832900393199363,      -6.547903238509202E-02,
        -0.3744576804916627,         -0.4417581235917372,         -0.8152491678999195,
        -0.7102685041338660,          0.6493006263927156,         -0.2718958414567380,
        -0.1659485415072040,         -0.7757475365189580,         -0.6088323588282737,
         0.8152491678999195,          0.3744576804916627,         -0.4417581235917372,
         0.3358108236422033,          0.8807282002850116,          0.3339894129272207,
        -0.8807282002850116,          0.3339894129272207,         -0.3358108236422033,
        -0.7757475365189580,          0.6088323588282737,          0.1659485415072040,
         0.6493006263927156,         -0.2718958414567380,         -0.7102685041338660,
         0.3358108236422033,         -0.8807282002850116,         -0.3339894129272207,
         0.8807282002850116,          0.3339894129272207,          0.3358108236422033,
         0.3339894129272207,          0.3358108236422033,          0.8807282002850116,
         0.2718958414567380,          0.7102685041338660,          0.6493006263927156,
         0.8152491678999195,         -0.3744576804916627,          0.4417581235917372,
        -0.6088323588282737,          0.1659485415072040,          0.7757475365189580,
     -6.547903238509202E-02,          0.1698622821349993,          0.9832900393199363,
        -0.1698622821349993,          0.9832900393199363,       6.547903238509202E-02,
         0.7102685041338660,          0.6493006263927156,          0.2718958414567380,
         0.9832900393199363,      -6.547903238509202E-02,          0.1698622821349993,
        -0.6493006263927156,         -0.2718958414567380,          0.7102685041338660,
         0.7757475365189580,          0.6088323588282737,         -0.1659485415072040,
        -0.4417581235917372,          0.8152491678999195,          0.3744576804916627,
         0.6088323588282737,          0.1659485415072040,         -0.7757475365189580,
        -0.8152491678999195,         -0.3744576804916627,         -0.4417581235917372,
        -0.1659485415072040,          0.7757475365189580,          0.6088323588282737,
         0.1698622821349993,         -0.9832900393199363,       6.547903238509202E-02,
        -0.6493006263927156,          0.2718958414567380,         -0.7102685041338660,
         0.9832900393199363,       6.547903238509202E-02,         -0.1698622821349993,
        -0.2718958414567380,          0.7102685041338660,         -0.6493006263927156,
      6.547903238509202E-02,          0.1698622821349993,         -0.9832900393199363,
        -0.3744576804916627,          0.4417581235917372,          0.8152491678999195,
        -0.7102685041338660,         -0.6493006263927156,          0.2718958414567380,
        -0.4417581235917372,         -0.8152491678999195,         -0.3744576804916627,
         0.7757475365189580,         -0.6088323588282737,          0.1659485415072040,
         0.1659485415072040,         -0.7757475365189580,          0.6088323588282737,
         0.3744576804916627,         -0.4417581235917372,          0.8152491678999195,
        -0.1996445654983616,         -0.8453065472879562,          0.4955793464008411,
        -0.8453065472879562,         -0.4955793464008411,          0.1996445654983616,
         0.3743841254383501,          0.7373115600658436,         -0.5623239191187172,
        -0.9933505989863832,         -0.1079949872221124,       3.989574198193997E-02,
        -0.1996445654983616,          0.8453065472879562,         -0.4955793464008411,
         0.8453065472879562,         -0.4955793464008411,         -0.1996445654983616,
         0.7619684846170787,         -0.4977712525855422,         -0.4142798674202900,
      3.989574198193997E-02,         -0.9933505989863832,         -0.1079949872221124,
         0.1079949872221124,       3.989574198193997E-02,          0.9933505989863832,
        -0.4142798674202900,         -0.7619684846170787,          0.4977712525855422,
        -0.7373115600658436,          0.5623239191187172,          0.3743841254383501,
         0.4310266798676661,         -0.2395403074803016,          0.8699634718391911,
         0.5623239191187172,          0.3743841254383501,         -0.7373115600658436,
         0.2395403074803016,          0.8699634718391911,         -0.4310266798676661,
         0.4977712525855422,          0.4142798674202900,          0.7619684846170787,
        -0.8699634718391911,         -0.4310266798676661,         -0.2395403074803016,
        -0.9933505989863832,          0.1079949872221124,      -3.989574198193997E-02,
        -0.4955793464008411,          0.1996445654983616,         -0.8453065472879562,
        -0.4955793464008411,         -0.1996445654983616,          0.8453065472879562,
         0.3743841254383501,         -0.7373115600658436,          0.5623239191187172,
         0.4955793464008411,         -0.1996445654983616,         -0.8453065472879562,
        -0.1079949872221124,       3.989574198193997E-02,         -0.9933505989863832,
         0.4142798674202900,         -0.7619684846170787,         -0.4977712525855422,
        -0.4310266798676661,         -0.2395403074803016,         -0.8699634718391911,
     -3.989574198193997E-02,          0.9933505989863832,         -0.1079949872221124,
        -0.5623239191187172,         -0.3743841254383501,         -0.7373115600658436,
        -0.7619684846170787,          0.4977712525855422,         -0.4142798674202900,
        -0.2395403074803016,         -0.8699634718391911,         -0.4310266798676661,
         0.7373115600658436,          0.5623239191187172,         -0.3743841254383501,
         0.1996445654983616,          0.8453065472879562,          0.4955793464008411,
        -0.8453065472879562,          0.4955793464008411,         -0.1996445654983616,
        -0.8699634718391911,          0.4310266798676661,          0.2395403074803016,
         0.4977712525855422,         -0.4142798674202900,         -0.7619684846170787,
         0.1996445654983616,         -0.8453065472879562,         -0.4955793464008411,
         0.8453065472879562,          0.4955793464008411,          0.1996445654983616,
         0.4955793464008411,          0.1996445654983616,          0.8453065472879562,
         0.4142798674202900,          0.7619684846170787,          0.4977712525855422,
         0.7373115600658436,         -0.5623239191187172,          0.3743841254383501,
        -0.4310266798676661,          0.2395403074803016,          0.8699634718391911,
        -0.1079949872221124,      -3.989574198193997E-02,          0.9933505989863832,
      3.989574198193997E-02,          0.9933505989863832,          0.1079949872221124,
         0.7619684846170787,          0.4977712525855422,          0.4142798674202900,
         0.9933505989863832,         -0.1079949872221124,      -3.989574198193997E-02,
        -0.4977712525855422,         -0.4142798674202900,          0.7619684846170787,
         0.8699634718391911,          0.4310266798676661,         -0.2395403074803016,
        -0.3743841254383501,          0.7373115600658436,          0.5623239191187172,
         0.4310266798676661,          0.2395403074803016,         -0.8699634718391911,
        -0.7373115600658436,         -0.5623239191187172,         -0.3743841254383501,
        -0.2395403074803016,          0.8699634718391911,          0.4310266798676661,
     -3.989574198193997E-02,         -0.9933505989863832,          0.1079949872221124,
        -0.4977712525855422,          0.4142798674202900,         -0.7619684846170787,
         0.9933505989863832,          0.1079949872221124,       3.989574198193997E-02,
        -0.4142798674202900,          0.7619684846170787,         -0.4977712525855422,
         0.1079949872221124,      -3.989574198193997E-02,         -0.9933505989863832,
        -0.5623239191187172,          0.3743841254383501,          0.7373115600658436,
        -0.7619684846170787,         -0.4977712525855422,          0.4142798674202900,
        -0.3743841254383501,         -0.7373115600658436,         -0.5623239191187172,
         0.8699634718391911,         -0.4310266798676661,          0.2395403074803016,
         0.2395403074803016,         -0.8699634718391911,          0.4310266798676661,
         0.5623239191187172,         -0.3743841254383501,          0.7373115600658436,
     -5.387602213028118E-02,         -0.7695747905192749,          0.6362798252629354,
        -0.7695747905192749,         -0.6362798252629354,       5.387602213028118E-02,
         0.2878105649135883,          0.6249952919365247,         -0.7256348694634922,
        -0.9573876030316789,         -0.1445794985827502,          0.2500114920374999,
     -5.387602213028118E-02,          0.7695747905192749,         -0.6362798252629354,
         0.7695747905192749,         -0.6362798252629354,      -5.387602213028118E-02,
         0.7795108915937734,         -0.3211077777687436,         -0.5378220569510883,
         0.2500114920374999,         -0.9573876030316789,         -0.1445794985827502,
         0.1445794985827502,          0.2500114920374999,          0.9573876030316789,
        -0.5378220569510883,         -0.7795108915937734,          0.3211077777687436,
        -0.6249952919365247,          0.7256348694634922,          0.2878105649135883,
         0.2317527335681867,         -0.3038875141677811,          0.9240903901765236,
         0.7256348694634922,          0.2878105649135883,         -0.6249952919365247,
         0.3038875141677811,          0.9240903901765236,         -0.2317527335681867,
         0.3211077777687436,          0.5378220569510883,          0.7795108915937734,
        -0.9240903901765236,         -0.2317527335681867,         -0.3038875141677811,
        -0.9573876030316789,          0.1445794985827502,         -0.2500114920374999,
        -0.6362798252629354,       5.387602213028118E-02,         -0.7695747905192749,
        -0.6362798252629354,      -5.387602213028118E-02,          0.7695747905192749,
         0.2878105649135883,         -0.6249952919365247,          0.7256348694634922,
         0.6362798252629354,      -5.387602213028118E-02,         -0.7695747905192749,
        -0.1445794985827502,          0.2500114920374999,         -0.9573876030316789,
         0.5378220569510883,         -0.7795108915937734,         -0.3211077777687436,
        -0.2317527335681867,         -0.3038875141677811,         -0.9240903901765236,
        -0.2500114920374999,          0.9573876030316789,         -0.1445794985827502,
        -0.7256348694634922,         -0.2878105649135883,         -0.6249952919365247,
        -0.7795108915937734,          0.3211077777687436,         -0.5378220569510883,
        -0.3038875141677811,         -0.9240903901765236,         -0.2317527335681867,
         0.6249952919365247,          0.7256348694634922,         -0.2878105649135883,
      5.387602213028118E-02,          0.7695747905192749,          0.6362798252629354,
        -0.7695747905192749,          0.6362798252629354,      -5.387602213028118E-02,
        -0.9240903901765236,          0.2317527335681867,          0.3038875141677811,
         0.3211077777687436,         -0.5378220569510883,         -0.7795108915937734,
      5.387602213028118E-02,         -0.7695747905192749,         -0.6362798252629354,
         0.7695747905192749,          0.6362798252629354,       5.387602213028118E-02,
         0.6362798252629354,       5.387602213028118E-02,          0.7695747905192749,
         0.5378220569510883,          0.7795108915937734,          0.3211077777687436,
         0.6249952919365247,         -0.7256348694634922,          0.2878105649135883,
        -0.2317527335681867,          0.3038875141677811,          0.9240903901765236,
        -0.1445794985827502,         -0.2500114920374999,          0.9573876030316789,
         0.2500114920374999,          0.9573876030316789,          0.1445794985827502,
         0.7795108915937734,          0.3211077777687436,          0.5378220569510883,
         0.9573876030316789,         -0.1445794985827502,         -0.2500114920374999,
        -0.3211077777687436,         -0.5378220569510883,          0.7795108915937734,
         0.9240903901765236,          0.2317527335681867,         -0.3038875141677811,
        -0.2878105649135883,          0.6249952919365247,          0.7256348694634922,
         0.2317527335681867,          0.3038875141677811,         -0.9240903901765236,
        -0.6249952919365247,         -0.7256348694634922,         -0.2878105649135883,
        -0.3038875141677811,          0.9240903901765236,          0.2317527335681867,
        -0.2500114920374999,         -0.9573876030316789,          0.1445794985827502,
        -0.3211077777687436,          0.5378220569510883,         -0.7795108915937734,
         0.9573876030316789,          0.1445794985827502,          0.2500114920374999,
        -0.5378220569510883,          0.7795108915937734,         -0.3211077777687436,
         0.1445794985827502,         -0.2500114920374999,         -0.9573876030316789,
        -0.7256348694634922,          0.2878105649135883,          0.6249952919365247,
        -0.7795108915937734,         -0.3211077777687436,          0.5378220569510883,
        -0.2878105649135883,         -0.6249952919365247,         -0.7256348694634922,
         0.9240903901765236,         -0.2317527335681867,          0.3038875141677811,
         0.3038875141677811,         -0.9240903901765236,          0.2317527335681867,
         0.7256348694634922,         -0.2878105649135883,          0.6249952919365247
};


static constexpr std::array<T,312> weights = {      2.678894659202470E-02,
      2.678894659202470E-02,
      2.678894659202470E-02,
      2.678894659202470E-02,
      2.678894659202470E-02,
      2.678894659202470E-02,
      2.678894659202470E-02,
      2.678894659202470E-02,
      2.678894659202470E-02,
      2.678894659202470E-02,
      2.678894659202470E-02,
      2.678894659202470E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      4.269044939387851E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      3.667266847936795E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.133414202330638E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.245640081450885E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02,
      4.092806020985292E-02
};
};
}}